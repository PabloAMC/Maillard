# Mahadevan & Farmer 2006 — EXTRACTION (file stem `mahadevan2006`)
J. Agric. Food Chem. **54**(19):7242-7250, DOI 10.1021/jf061102x (printed p. 7242). Crossref 2026-11-03: "Key Odor
Impact Compounds in Three Yeast Extract Pastes", authors Mahadevan, Farmer; JAFC 54, pp. 7242-7250, published online
2006-08-17. **Matches the PDF and the brief.** No naming caveat.

### Three commercial *S. cerevisiae* yeast-extract pastes (types A [batches X, Y], B, C; two suppliers; factory-made), SDE (Likens-Nickerson, 2 h) and dynamic headspace (60 °C, 30 min, Tenax) → GC-O with **odour scores on a 1-5 scale** (undiluted and 1/10), GC-MS / GC-IR identification; precursors (thiamine, 13 ribonucleotides / nucleosides, ribose, xylose, fructose, glucose) by HPLC in **mg/100 g dry weight**. **No odorant concentration in any unit. No heating series, no time, no kinetics; factory processing conditions only sketched (70 °C under vacuum).**

**Read 2026-11-03 from `data/articles/mahadevan2006.pdf`** (9 pp., ACS typeset, text layer present). Prose from
`pdftotext -layout` (scratchpad `ingest/mahadevan2006_work/text.txt`). **Table 1 (p. 7245) and Table 2 (p. 7246)
read by eye from 150 dpi renders** (`p-4.png`, `p-5.png`); **Tables 3-5 (p. 7247) from a 200 dpi render**
(`p-6.png`). The text layer agrees with the renders cell for cell for Tables 1, 2 and 5. No figures in the paper.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Food Science, Queen's University Belfast (Mahadevan, PhD student) and Agriculture, Food
and Environmental Science Division, **Agri-Food and Biosciences Institute (AFBI)**, Newforge Lane, Belfast;
corresponding author **Linda Farmer** (linda.farmer@afbini.gov.uk). **Existing id `afbi_belfast_farmer`** ("AFBI
Belfast (Farmer)") **fits** (P ~0.95): same corresponding author, institution and address as the key's evidence
(Aliani 2013); the nucleotide and sugar HPLC methods are Aliani & Farmer 2002 / 2005 (refs 19, 20) [C]. Same-lab
dossier: `aliani2013_extraction.md`. `farmer1989/1990/1992_extraction.md` are Farmer's AFRC Bristol (Mottram-era)
work, which the key leaves unmapped; same person, different institution. **Independent of DFA Garching**
(`munch1997`, `munch1998`) and of BTBU (`alim2018`, `alim2019`).

## 1. What this paper answers for the model

| model need | what Mahadevan 2006 gives | tag |
|---|---|---|
| MFT, FFT, methional amounts in yeast extract | **none**: GC-O odour scores (1-5, two common assessors) only | [M] / [!] |
| which odorants matter in YE | Table 1 (SDE, 31 odours, 3 pastes), Table 2 (headspace, 31 odours, 2 batches of A) | [M] |
| sulfur odorants present | MFT, FFT, 2-methyl-3-(methyldithio)furan, 2-methyl-3-(methyltrithio)furan, bis(2-methyl-3-furyl) disulfide, DMTS, DMTetraS, 2-methylthiophene, 2-methyltetrahydrothiophen-3-one, trithiolane, dithiolanone | [M] (scores) |
| Strecker aldehydes | methional, phenylacetaldehyde (scores); 2-/3-methylbutanal **not among odours** | [M] |
| pyrazines | 2-methyl-5-isopropyl- (tentative), trimethyl-, 2-ethyl-3,6-dimethyl-, 2-ethyl-3,5-dimethylpyrazine (scores) | [M] |
| furanones (HDMF, sotolon) | **not reported** | — |
| precursor composition | Table 5: thiamine, 5'-IMP/GMP/CMP/UMP/AMP/ADP/ATP, 2'/3'-AMP, 3'-GMP, nucleosides, ribose, xylose, fructose, glucose; mg/100 g DW, duplicate, recoveries | [M] |
| free amino acids, cysteine, peptides | **not measured** | — |
| heating conditions of the odorants | factory: autolysate → liquor at 70 °C under vacuum, evaporation to paste; B evaporated longer; C unknown; SDE adds 2 h at ~100 °C | [M] / [Q] |
| pH | "pH of 5.1" for all three diluted pastes (dilution not stated) | [M] |

**Answers.**
1. **Which odorants?** (confidence ~95 % on the list). Key odours in all three pastes: MFT (LRI 1310), FFT (1432),
   methional (1455), β-damascenone (1821), 2-methoxyphenol (1861), plus 2-methyl-3-(methyldithio)furan (1669), DMTS
   (1373), 1-octen-3-one (1298, by LRI + odour only, not seen by MS). New as YE key odorants per the authors [Q]:
   2-methylthiophene, DMTS, 1-octen-3-one, 1-methylthio-3-pentanone, 5-methylfurfural, the methyl di- and trithio
   furans. **No HDMF, sotolon, 2-/3-methylbutanal, 2,3-butanedione** among the odours (contrast `munch1998`, where
   these are SIDA-quantified with OAV 99-3640) [!].
2. **Units?** (confidence ~99 %). Odorants: **GC-O intensity scores 0-5** (mean of two assessors, duplicate) at
   undiluted and 1/10 dilution — not FD factors, not concentrations, not relative peak areas. Bromobenzene was added as
   an SDE internal standard but **no IS-based amount is printed** in the paper; the full volatile list is in the SI
   (sec. 8). Precursors: **absolute, mg/100 g dry weight**, HPLC with external standards and spike recoveries,
   duplicate, "without correction for recovery" [V].
3. **Heating conditions?** (confidence ~90 % that nothing simulable is given). The pastes are factory products:
   autolysis, then liquor (~45 % solids) by heating "under vacuum at 70 °C" [V], then evaporation to 75-83 % solids;
   B "thermally evaporated for a longer period" [V]; no times. Type C (other supplier): nothing. The SDE extracts carry
   a further 2 h of boiling; headspace collection 30 min at 60 °C.
4. **Non-volatiles?** Thiamine 1.12-3.00, 5'-IMP 88-496, 5'-GMP 98-138, ribose 6.0-13.2, xylose 1.1-4.2, fructose
   10.9-21.2, glucose 1.6-2.4 mg/100 g DW (Table 5) [M]. **No amino acids, Cys, peptides.**
5. **A second lab's quantitative YE odorant data comparable with Münch 1998?** **No** (confidence ~97 %). Same
   odorants (MFT, FFT, methional, 2-AP, phenylacetaldehyde) but no amounts, and the only route to amounts (the SI
   volatile table) is at best bromobenzene-relative peak areas (P ~0.7 it carries no absolute concentrations; SI not on
   disk). What it does give: **directional between-paste contrasts** (B > A for MFT, DMTS, methyldithiofuran; C has
   ×5 the 5'-IMP of A/B and the strongest MFT score) and a **second lab's YE precursor composition** (nucleotides,
   sugars, thiamine), which `munch1997/1998` lack.
6. **Mechanism claims** are [P]: IMP → ribose → MFT / methyldithiofuran (cites Mottram & Madruga [C]); di-/trisulfides
   from methanethiol oxidation; thiazoles from Maillard rather than thiamine because C has least thiamine.

## 2. Conditions

- **Material** [M]: three *S. cerevisiae* YE pastes for savoury foods. A (batch X 75.0 %, batch Y 76.9 % total solids)
  and B (82.7 %) from the same supplier and liquor; C (80.6 %) from another supplier. "Pure": no added spices [Q].
- **pH** [M]: 5.1 for all three "when diluted" (dilution and T not stated). No buffer.
- **Process** [M]: autolysate heated under vacuum at 70 °C to liquor (~45 % solids), evaporated (same conditions) to
  paste; B evaporated longer. No times, no temperatures beyond 70 °C.
- **SDE** [M]: 60 g paste in 200 mL water, pentane/ether 1:1 (20 mL) with bromobenzene IS (20 µL of 1000 ng/mL), 2 h;
  duplicate extractions; A (batch Y), B, C.
- **Dynamic headspace** [M]: GC-MS: 40 g paste in 100 mL water, 60 °C, 30 min, N2 50 mL/min onto 2.6 mg (sic) Tenax
  GC; GC-O: **0.5 g dry weight in 100 mL** (A batches X and Y). Duplicate collections (A-Y, B).
- **Replicates** [M]: GC-O 4 assessors (2 for C); scores analysed by REML for the two common assessors.

## 3. Analytical methods

- **GC-MS** HP 6890/5973, CPWax52CB 50 m and CPSil8CB 50 m; GC-IR (HP 5965B) for confirmation. Identity: MS + LRI +
  odour of authentic standard ("MS, LRI, O") or literature only (lower case, tentative) [M].
- **GC-O** [M]: undiluted and 1/10; scores 1 (very weak) to 5 (very strong), 60 min; methods say only odours with
  mean ≥ 0.25 (diluted) and ≥ 1.2 (undiluted) are reported, but the table footer counts odours "≥ 0.05" [!].
- **Precursors** [M]: thiamine as thiochrome by HPLC (recovery 102 ± 4 %); ribonucleotides by perchloric acid
  extraction + ion-pair HPLC, two programs (IMP / GMP co-elute in program 1); sugars by HPLC with post-column
  derivatisation, cellobiose IS. Duplicates. ± in Table 5 = SEM of the ANOVA; recoveries ± SD (duplicate).

## 4. Data

### 4a. Table 1 (p. 7245): SDE GC-O mean odour scores, pastes A (batch Y) / B / C [M]
Undiluted | diluted 1/10; P for paste effect (two common assessors); identification method. Rows with an identified
compound; eight rows assigned "unknown" (LRI 1391 partly, 1620, 1693, 1795, 2394, and co-eluting unknowns) omitted
except LRI 1620 / 1693, which the authors discuss.

| LRI | compound | A / B / C undiluted | A / B / C 1/10 | P | ident. |
|---|---|---|---|---|---|
| 1101 | 2-methylthiophene | 2.9 / 1.8 / 2.8 | 0.1 / 0.5 / 0.0 | ns | MS, LRI, O |
| 1232 | 2-pentylfuran | 0.0 / 0.6 / 0.8 | 0.0 / 0.0 / 1.0 | ns | MS, LRI, O |
| 1285 | 3-hydroxy-2-butanone + bis(methylthio)methane | 0.0 / 2.3 / 1.3 | 0.0 / 1.0 / 0.5 | <0.001 | MS, LRI, O |
| 1298 | 1-octen-3-one | 0.0 / 3.1 / 3.3 | 0.0 / 2.5 / 3.1 | <0.001 | LRI, O |
| 1310 | **2-methyl-3-furanthiol** | 2.1 / 3.5 / 3.4 | 0.6 / 2.0 / 2.8 | <0.001 | MS, LRI, O |
| 1335 | 2-acetyl-1-pyrroline (tent.) | 2.6 / 2.9 / 2.5 | 0.6 / 1.4 / 1.0 | ns | lri, o |
| 1356 | 2-methyl-3-(methylthio)furan (tent.) | 0.8 / 0.8 / 2.2 | 0.4 / 0.0 / 0.0 | ns | ms, lri, o |
| 1373 | dimethyl trisulfide | 0.4 / 3.3 / 2.1 | 0.0 / 1.1 / 0.9 | <0.001 | MS, LRI, O |
| 1391 | 2-nonanone + unknown | 0.0 / 1.4 / 2.2 | 0.0 / 0.0 / 1.0 | 0.008 | MS, LRI, O |
| 1415 | 2-methyl-5-isopropylpyrazine (tent.) | 3.1 / 2.5 / 2.3 | 0.0 / 0.0 / 0.6 | ns | ms, lri |
| 1432 | **2-furanmethanethiol (FFT)** | 3.8 / 3.6 / 4.2 | 3.3 / 1.5 / 2.4 | 0.063 | MS, LRI, O |
| 1455 | **methional** | 3.1 / 3.8 / 3.8 | 1.5 / 2.3 / 2.5 | 0.014 | MS, LRI, O |
| 1468 | furfural + unknown | 0.0 / 0.5 / 1.7 | 0.0 / 0.3 / 0.0 | ns | MS, LRI, O |
| 1523 | 2-methyltetrahydrothiophen-3-one + unknown | 0.0 / 1.3 / 3.8 | 0.0 / 0.4 / 3.0 | <0.001 | MS, LRI, O |
| 1570 | 5-methylfurfural + 1-methylthio-3-pentanone | 1.8 / 1.0 / 2.6 | 0.5 / 0.0 / 0.3 | ns | MS, LRI, O |
| 1586 | cis-3,5-dimethyl-1,2,4-trithiolane (tent.) | 0.0 / 1.9 / 0.6 | 0.0 / 0.0 / 0.0 | <0.001 | MS, LRI |
| 1620 | unknown (potato, savoury) | 3.3 / 2.5 / 2.5 | 2.1 / 0.5 / 0.8 | <0.001 | — |
| 1644 | phenylacetaldehyde | 3.0 / 3.0 / 1.3 | 0.3 / 0.0 / 0.0 | ns | MS, LRI, O |
| 1669 | **2-methyl-3-(methyldithio)furan** | 0.8 / 2.3 / 2.8 | 0.8 / 1.8 / 2.4 | <0.001 | MS, LRI, O |
| 1693 | unknown (potato, roasted) | 2.9 / 2.8 / 4.0 | 1.5 / 1.4 / 1.3 | ns | — |
| 1735 | dimethyl tetrasulfide | 0.0 / 0.9 / 1.3 | 0.0 / 1.0 / 0.5 | 0.018 | LRI, O |
| 1755 | 2-acetyl-2-thiazoline | 0.0 / 0.6 / 1.2 | 0.0 / 0.0 / 0.4 | ns | LRI, O |
| 1782 | 5-methyl-2-thiophenecarboxaldehyde + unknown | 1.8 / 0.6 / 1.1 | 0.0 / 0.0 / 0.0 | ns | MS, LRI, O |
| 1821 | trans-β-damascenone | 2.9 / 2.6 / 3.2 | 1.3 / 1.3 / 2.4 | ns | MS, LRI, O |
| 1861 | 2-methoxyphenol | 2.8 / 2.1 / 3.2 | 1.1 / 0.9 / 1.0 | ns | MS, LRI, O |
| 1932 | 2-phenylbut-2-enal | 0.0 / 0.8 / 0.7 | 0.0 / 0.6 / 0.5 | 0.09 | MS, LRI, O |
| 1976 | 2-methyl-3-(methyltrithio)furan | 1.5 / 2.3 / 1.8 | 0.0 / 0.6 / 0.8 | ns | MS, LRI, O |
| 2022 | 2-methylphenol | 0.5 / 0.5 / 1.0 | 0.0 / 0.6 / 0.8 | ns | MS, LRI, O |
| 2145 | bis(2-methyl-3-furyl) disulfide | 0.4 / 0.8 / 0.9 | 0.4 / 1.0 / 0.8 | ns | MS, LRI, O |
| — | overall mean score | 1.3 / 1.8 / 2.2 | 0.5 / 0.8 / 1.1 | | |
| — | no. odours detected | 19 / 31 / 31 | 14 / 21 / 25 | | |

### 4b. Table 2 (p. 7246): dynamic-headspace GC-O mean scores, paste A batch X / batch Y (4 assessors) [M]
Sulfur, Strecker and pyrazine rows (X / Y; P if < 0.10): 1310 MFT 1.5 / 3.0 (0.069); 1338 2-methyl-3-(methylthio)furan
(tent.) 0.0 / 1.0; 1375 DMTS 3.4 / 3.1; 1400 trimethylpyrazine + unknown 1.3 / 1.6; 1437 FFT 0.0 / 1.9 (0.022);
1443 2-ethyl-3,6-dimethylpyrazine 1.3 / 1.3; 1456 methional + furfural + 2-ethyl-3,5-dimethylpyrazine 3.1 / 3.0;
1485 2-furfuryl methyl sulfide + a pyrazine 0.4 / 1.5; 1525 2-methyltetrahydrothiophen-3-one 1.2 / 1.3; 1605
trans-3,5-dimethyl-1,2-dithiolan-4-one 0.0 / 1.0; 1642 phenylacetaldehyde 2.3 / 2.0; 1648 3-methylthiothiophene (tent.)
1.6 / 2.0; 1667 2-methyl-3-(methyldithio)furan 2.8 / 2.8; 1769 2-acetylthiophene (tent.) 0.0 / 1.5; 1937
2-thiophenemethanol + benzothiazole 0.1 / 2.3 (0.021); 1979 2-methyl-3-(methyltrithio)furan 0.0 / 2.8 (0.022).
Other rows: 956 2-ethylfuran (tent.) 2.5 / 1.3; 1049 2,3-pentanedione 0.0 / 1.5; 1297 1-octen-3-one 2.1 / 1.8
(0.008); 1340 2-AP (tent.) 0.7 / 1.5; 1810 (E,E)-2,4-decadienal 0.0 / 1.3; 1824 β-damascenone (tent.) 1.4 / 1.3; 1861
2-methoxyphenol 0.0 / 2.5 (0.030). Overall mean 1.0 / 1.7; odours detected 18 / 31. GC-MS was not run on batch X
(insufficient sample): batch X identities are borrowed from batch Y [M] [!].

### 4c. Tables 3 / 4: sensory descriptors only (beef, root vegetables, yeast extract...; C > A, B; Y > X) [Q].

### 4d. Table 5 (p. 7247): precursors, mg/100 g dry weight, mean of duplicates, not recovery-corrected [M]
Letters: Fisher LSD groups (p < 0.05); sig. *** P ≤ 0.001, ** ≤ 0.01, * ≤ 0.05; SEM; recovery % ± SD.

| analyte | paste A | paste B | paste C | sig | SEM | recovery % |
|---|---|---|---|---|---|---|
| thiamine | 2.58 b | 3.00 b | 1.12 a | *** | 0.16 | 102 ± 4 |
| 5'-IMP | 93 a | 88 a | 496 b | *** | 4.63 | 98 ± 2 |
| inosine | 71 b | 74 b | 27 a | *** | 0.39 | |
| 5'-GMP | 99 a | 98 a | 138 b | *** | 1.24 | |
| guanosine | 213 b | 291 c | 124 a | *** | 1.10 | 103 ± 3 |
| 5'-CMP | 226 a | 340 a | 480 b | *** | 7.20 | |
| 5'-UMP | 230 b | 203 b | 93 a | *** | 1.89 | 87 ± 2 |
| 5'-ATP | 68 | 72 | 65 | NS | 5.63 | 100 ± 2 |
| 5'-ADP | 436 a | 585 b | 501 a | ** | 9.71 | 66 ± 8 |
| 5'-AMP | 150 b | 138 b | 80 a | *** | 1.28 | |
| adenosine | 469 b | 625 c | 321 a | *** | 6.32 | |
| 3'-AMP | 344 b | 361 b | 206 a | ** | 7.69 | |
| 2'-AMP | 1566 a | 1870 b | 1528 a | *** | 10.97 | |
| 3'-GMP | 686 a | 831 b | 697 a | ** | 6.76 | |
| ribose | 11.08 b | 6.04 a | 13.17 b | * | 1.02 | 72 ± 0.3 |
| xylose | 4.17 b | 1.07 a | 1.33 a | ** | 0.17 | 98 ± 1 |
| fructose | 10.92 a | 16.45 ab | 21.15 b | ** | 0.89 | 97 ± 2 |
| glucose | 1.92 a | 1.62 a | 2.40 b | ** | 0.06 | 87 ± 3 |

Derived [Z]: total ribonucleotides + nucleosides (all 13 rows) A 4651, B 5576, C 4756 mg/100 g DW (≈ 5 % of dry
matter, a large potential ribose pool next to 6-13 mg/100 g free ribose); 2'+3' isomers / that total A 0.558, B 0.549, C 0.511 vs printed 56 / 61 / 51 %: **A and C reproduce
(only with nucleosides in the denominator), B does not (55 vs 61 %)** [!]. Free reducing sugars total
25-38 mg/100 g DW, i.e. ~0.03 % of DW: **free sugar is scarce; nucleotide-bound ribose dominates** [Z].

## 5. Flags and internal consistency

1. [!] **No odorant amounts** anywhere in the main text; GC-O scores are ordinal intensities of two assessors.
2. [!] Reporting threshold: methods say mean ≥ 0.25 (diluted) / ≥ 1.2 (undiluted), footer counts "≥ 0.05", and many
   listed rows are 0.0 for the paste in question (they qualified in another paste).
3. [!] Table 5 does not say which batch "paste A" is (X or Y); the SDE and sensory work used batch Y (P ~0.7 it is Y).
4. [!] 5'-CMP: A 226 and B 340 share letter "a" with SEM 7.20 (difference ~16 SEM); a lettering slip or SEM is not the
   per-mean error. Likewise ribose A 11.08 vs C 13.17 share "b" with SEM 1.02 — plausible; CMP is not.
5. [!] Ribose recovery 72 ± 0.3 % and ADP 66 ± 8 %, not corrected: true ribose ~1.4 × the printed value [Z].
6. [!] Headspace GC-O used 0.5 g DW / 100 mL, GC-MS 40 g paste / 100 mL (~60 × more); batch X identities inferred from Y.
7. [!] B: 2'+3' share of total ribonucleotides 55 % from Table 5 vs 61 % in the text [Z]. ("2.6 mg" Tenax is implausibly small; immaterial here.)
8. pH 5.1 for all pastes is printed once, without dilution; consistent with the text's later "pH of all three pastes was
   the same".
9. Text ↔ table: the P values quoted in the text (P < 0.001 for B > A odours, LRI 1620 A > B) match Table 1.

## 6. Existing records in the repo that cite this paper

grep 2026-11-03 for `jf061102x`, `mahadevan`, "Three Yeast Extract Pastes", "Farmer ... yeast" across `data/`,
`docs/`, `src/`, `scripts/`, `tasks/`:
- `data/keys/pdf_dossier_aliases.yml` l. 98 (pending queue line): authors, year, JAFC 54:7242, DOI, "key odour-impact
  compounds in three yeast-extract pastes" — **match**; "a second lab for yeast-extract odorants" — **true as to lab,
  not as to quantitative odorants** (scores only).
- `data/research_corpus/raw/Gemini_deep_research_04.md` l. 252: bibliography item 44 (title + DOI) — **match**; item 44
  is not cited in the report body (grep), so no claim to check.
- No ledger, registry, benchmark, prior or code record cites it.

## 7. Proposed ledger entries (not written into the ledgers)

```yaml
# sulfur.yml
  # ---- mahadevan2006 : Mahadevan & Farmer 2006 (GC-O of three commercial yeast-extract pastes) (ingested 2026-11-03) ----
  - id: mahadevan2006_tab1_tab2_gco_scores_yeast_extract_pastes
    dossier: mahadevan2006_extraction.md
    paper: "Mahadevan & Farmer 2006 J Agric Food Chem 54:7242"
    lab: "AFBI Belfast (Farmer)"
    system: "three commercial S. cerevisiae yeast-extract pastes (A batches X/Y, B same supplier longer evaporation, C other supplier), 75-83 % solids"
    medium: "factory-made pastes; SDE 60 g in 200 mL water 2 h; dynamic headspace 40 g (MS) or 0.5 g DW (GC-O) in 100 mL water, 60 C 30 min"
    temperature_c: []
    ph: {initial: 5.1, measured_during: false, final: null}
    time_min: []
    species: ["MFT", "FFT", "2-methyl-3-(methyldithio)furan", "2-methyl-3-(methyltrithio)furan", "bis(2-methyl-3-furyl) disulfide", "dimethyl trisulfide", "dimethyl tetrasulfide", "2-methylthiophene", "2-methyltetrahydrothiophen-3-one", "cis-3,5-dimethyl-1,2,4-trithiolane", "2-acetyl-2-thiazoline"]
    kind: ranking
    location: "Table 1 (SDE, pastes A/B/C, undiluted and 1/10) and Table 2 (headspace, batches X/Y) (dossier secs. 4a, 4b)"
    status: digitised_in_dossier
    quant_class: other
    simulable: false
    condition_quote: "Mean odor scores on a scale of 1 (very weak) to 5 (very strong)"
    flags: [matrix_food, condition_uncertain, no_time_temperature_history, ph_dilution_not_stated, "reporting_threshold_text_vs_footer", "batch_x_identities_from_batch_y", "also_relevant: aroma_other"]
    note: "GC-O intensity scores (two common assessors, duplicate), not FD factors or amounts. Directional: B (longer evaporation) > A for MFT, DMTS, methyldithiofuran; C (5x IMP) strongest MFT. Second-lab presence check for Muench 1998 odorants, not a quantitative comparison."
# aroma_other.yml
  - id: mahadevan2006_tab5_yeast_extract_precursors
    dossier: mahadevan2006_extraction.md
    paper: "Mahadevan & Farmer 2006 J Agric Food Chem 54:7242"
    lab: "AFBI Belfast (Farmer)"
    system: "three commercial S. cerevisiae yeast-extract pastes A, B, C (batch of A not stated)"
    medium: "pastes as sold, before any further heating"
    temperature_c: []
    ph: {initial: 5.1, measured_during: false, final: null}
    time_min: []
    species: [thiamine, "5'-IMP", "5'-GMP", "5'-CMP", "5'-UMP", "5'-AMP", "5'-ADP", "5'-ATP", "2'-AMP", "3'-AMP", "3'-GMP", inosine, guanosine, adenosine, ribose, xylose, fructose, glucose]
    kind: other
    location: "Table 5 (dossier sec. 4d)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "quantities are quoted without correction for recovery"
    flags: [matrix_food, dry_weight_basis, composition_incomplete, "ribose_recovery_72pct_uncorrected", "cmp_lettering_inconsistent_with_sem", "batch_of_paste_a_not_stated", "also_relevant: sulfur", "also_relevant: trunk"]
    note: "mg/100 g DW, HPLC external standards, duplicate, SEM. Free reducing sugars only ~0.03 % DW while ribonucleotides + nucleosides ~5 % DW: a nucleotide-ribose pool for MFT. No amino acids, Cys, peptides. Complements munch1997/1998 (amino acids, no sugars or nucleotides)."
```

No trunk.yml entry: no sugar, furan or browning data set beyond the composition row (flagged also_relevant).

## 8. Supplementary material

Yes: "Table of volatile compounds collected by simultaneous distillation and extraction from three types of yeast
extract paste samples", free at pubs.acs.org (p. 7249). It is the only place a per-compound amount could be
(bromobenzene-relative, at best). **Not on disk.** Worth fetching only if a relative-abundance YE volatile profile is
wanted (P ~0.3 it holds numbers beyond presence / LRI). The paper says follow-up papers on processing conditions and
precursors were "in preparation" [Q]; not identified here.
