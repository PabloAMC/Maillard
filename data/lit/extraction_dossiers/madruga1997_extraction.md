# Madruga 1997 — EXTRACTION (5'-IMP vs cysteine vs thiamine added to beef, 140 °C / 30 min: sulfur-substituted furans in the headspace)

**Source on disk:** `data/articles/madruga1997.pdf` (a 6-page scan, 150 ppi JPEG per page, **no text layer**,
so there was no `pdftotext` cross-check). Read and checked by eye on 2026-10-09 from the page images,
re-rendered at 300 and 600 dpi and zoomed row by row for Table 1. The body, summary and table are in
English; only the "Resumo" is in Portuguese. Column order in Table 1 is **Meat + 5'-IMP | Meat + THI |
Meat + CYS | Meat BLK** (not control-first); every value below was read against that header. Written for
the thiamine-vs-cysteine question raised by the composition sweep (thiamine predicted to add almost nothing
to MFT at 140 °C/5 min), not for the core fit.

| field | value |
|---|---|
| Title | "Studies on some precursors involved in meat flavour formation" (Portuguese title in the Resumo: "Estudos de alguns precursores envolvidos na formação do aroma cárneo") |
| Authors | M.S. Madruga (sole author; Departamento de Tecnologia Química e de Alimentos, DTQA/UFPB, João Pessoa, Paraíba) |
| Venue | Ciência e Tecnologia de Alimentos (Campinas) 1997, 17(2), 148-153 (running footer: "Ciênc. Tecnol. Aliment., 17 (2):148-153, mai-ago. 1997"); received 27/12/96, accepted 02/06/97 |
| DOI | 10.1590/S0101-20611997000200016 (not printed on the scan; from the request, consistent with volume, issue and pages) |

## 1. Methods

**Matrix and additions (§2.1, p. 148).** "Portions (100g) of minced beef *M. Psoas major*", from a local
supplier, chopped and blended "with 12 ml water containing different meat flavour precursors". One portion
each:

| treatment (Table 1 column) | addition per 100 g beef |
|---|---|
| Meat + 5'-IMP | 5'-IMP 2.7 g |
| Meat + THI | thiamine 1.0 mg |
| Meat + CYS | cysteine 39.3 mg |
| Meat BLK | no precursor (blank) |

Author's statement: "This resulted in an increase of approximately ten times in the concentration of these
precursors compared with normal meat." The native levels behind the "ten times" are not printed. The salt
forms of IMP and thiamine are not printed. No sugar was added in any arm (IMP is the pentose source: "the
principal source of pentose sugar in muscles", p. 148).

Molar doses, **derived here** (form assumptions stated): cysteine 39.3 mg ÷ 121.16 g/mol = 0.324 mmol per
100 g; thiamine 1.0 mg ÷ 265.35 g/mol (cation) to ÷ 337.27 g/mol (hydrochloride) = 0.0030 to 0.0038 mmol
per 100 g. Cysteine:thiamine ≈ 86 to 109 mol/mol. Implied native levels if "ten times" is taken literally
(derived): thiamine ≈ 0.1 mg/100 g, cysteine ≈ 3.9 mg/100 g; the thiamine figure is just above the
0.01–0.08 mg/100 g raw beef range in Lombardi-Boccia 2005.

**pH (§2.1).** "The addition of 5'-IMP resulted in a drop in pH of 1.1 pH units; therefore the pH was
adjusted, before heating, to 5.6 by adding 1 M sodium hydroxide." The pH of the THI, CYS and blank
portions is not printed (only the IMP portion is said to have been adjusted, to 5.6).

**Cooking (§2.1).** Left overnight in a refrigerator, then "heated in glass bottles in an autoclave at
140°C for 30min". Closed bottles, so a wet, sealed system. One temperature, one time.

**Volatile measurement (§2.2, p. 148-149).** Headspace from the cooked meat in a 250 ml flask at 60 °C with
agitation, swept by oxygen-free nitrogen 40 ml/min for 2 h onto Tenax-GC; GC-FID with odour port (DB-5,
30 m × 0.32 mm) and GC-MS (HP5988A). "Quantization of volatiles was based on peak area integration of the
GC-MS chromatograms using 1,2-diclorobenzene as internal standard"; 65 ng in 1 ml ethanol added to the trap
after collection (p. 149). Table 1 footnote 1: "concentrations obtained by comparing GC/MS peak areas with
the area of 65 ng dichlorobenzene added to Tenax trap as internal standard"; p. 150 adds that the areas
were in the total ion chromatogram. So: **one internal standard, TIC, response factor 1 for every
compound, no calibration**. The author: the method "only gives an indication of the approximate quantity
present in collected headspace volatiles". Units "ng/100g meat". Footnote 2: "tr, trace (.2 ng/100g
meat); nd, not detected; + present in significant amounts but quantitation confunded by large adjacent
peak".

**Replicates and statistics.** Table 1 caption: "Each value is the mean of three analysis and the standard
deviation is given in parentheses." §2.2.1: "Usually no less than three headspace collections were
performed for each sample." p. 150: "Values were averaged over at least three replicates for each
system." This reads as **repeat headspace collections from one cooked portion per treatment** (technical
replicates), not independent cooks; not stated outright. No significance test is reported anywhere,
although the Summary says sulfur furans were "significantly affected".

**GC-O.** Four assessors sniffed the GC effluent (§2.2.2) and four assessed whole-sample odour (§2.2.4).
Only descriptive results; no per-compound intensities or dilution factors are printed.

## 2. Findings that matter

### 2a. Table 1 (p. 150), S-substituted furans, ng/100 g meat, mean (SD)

| compound | LRI | Meat + 5'-IMP | Meat + THI | Meat + CYS | Meat BLK |
|---|---|---|---|---|---|
| 2-methyl-3-furanthiol (MFT) | 873 | 13 (5.9) | tr | 3 (1.0) | tr |
| 2-furylmethanethiol (2-furfurylthiol, FFT) | 913 | 250 (21.5) | 207(38) | 124(24) | 142(35) |
| 2-methyl-3-furyl methyl disulfide | 1175 | 23 (6.1) | 19 (4.9) | 13 (3.3) | tr |
| 2-furylmethyl methyl disulfide | 1220 | 74 (12.9) | 19 (3.6) | 23 (14.9) | 21 (15.7) |
| 2-methyl-3-furyl methyl trisulfide | 1392 | 3 (1.3) | 4 (0.2) | 2 (1.3) | tr |
| 2-furylmethyl methyl trisulfide | 1451 | 55 (12.2) | 19 (3.6) | 9 (1.4) | 30 (11.9) |
| 1-(2-methyl-3-furyldithio)-2-propanone | 1466 | tr | nd | tr | nd |
| bis(2-methyl-3-furyl) disulfide | 1535 | tr | tr | nd | nd |
| 2-(2-methyl-3-furyldithio)-3-pentanone | 1584 | tr | nd | nd | nd |
| 2-methyl-3-(2-furylmethyldithio)furan | 1635 | 8 (1.9) | 5 (0.5) | 2 (0.2) | 4 (1.4) |
| bis(2-furylmethyl) disulfide | 1687 | (blank cell) | 20 (12.5) | (blank cell) | 30 (5.9) |
| bis(2-furylmethyl) trisulfide | 1932 | 13 (8.2) | 6 (0.9) | 17 (1.6) | 7 (1.9) |

The bis(2-furylmethyl) disulfide row has empty cells under IMP and CYS; neither "nd", "tr" nor "+" is
printed there, so those two values are unknown.

Text (p. 151): bis(2-methyl-3-furyl) disulfide "was only detected in meat with added precursors"; the table
shows it at trace in IMP and THI only, nd in CYS as well as in the blank. MFT formation is attributed to
hydroxyfuranone and dicarbonyls "(from 5'-IMP)" reacting with "hydrogen sulfide (from cysteine or
thiamine)". The 2-furylmethyl compounds "were formed in all systems, including meat blank, and the amounts
formed were generally unaffected by the addition of precursors".

### 2b. Other sulfur rows from Table 1 (p. 149-150), for context

| compound | LRI | IMP | THI | CYS | BLK |
|---|---|---|---|---|---|
| 2-methylthiophene | 779 | 16 (6) | 13 (1.2) | 10 (6.1) | 14 (3.3) |
| 2-formylthiophene | 1002 | 191 (17.5) | 7 (3.5) | 602 (67) | 2 (0.7) |
| 2-acetylthiophene | 1087 | 15 (2.3) | nd | 59 (22) | nd |
| 3-methyl-1,2-dithiolan-4-one | 1072 | 45 (21.9) | 42 (3.8) | 83 (45.4) | 43 (9.4) |
| 2-acetylthiazole | 1020 | 227 (41) | 194 (23) | 269 (41) | 165 (36) |
| dimethyl disulfide | 763 | 13 (4.9) | 14 (5.1) | 22 (2.9) | 12 (1.7) |
| dimethyl trisulfide | 967 | 286 (19.2) | 269(27) | 310 (27.9) | 267 (40.8) |
| 3-(methylthio)propanal | 907 | 45 (10.4) | 9 (2.6) | 30 (2.9) | tr |
| 2-furfural | 828 | 36 (17.7) | 12 (7.2) | 18 (6.4) | 50 (26.5) |

### 2c. Derived comparisons (all **derived here** from Table 1)

- **MFT itself:** THI = tr (≈ 0.2 ng/100 g by footnote 2) = blank; CYS 3 (1.0); IMP 13 (5.9). Thiamine
  gave no MFT above the blank; cysteine gave a small amount; IMP the most.
- **Methyl 2-methyl-3-furyl disulfide (MFT + methanethiol product):** THI 19 (4.9) vs CYS 13 (3.3) vs
  blank tr. Thiamine ≥ cysteine. Welch t on n = 3 each, assuming the SDs are sample SDs:
  SE = √(4.9²/3 + 3.3²/3) = √11.6 = 3.41; t = (19 − 13)/3.41 = 1.76, about 3.5 df, p ≈ 0.16. Not
  distinguishable; both clearly above blank.
- **2-methyl-3-furyl pool** (MFT + methyl disulfide + methyl trisulfide, trace counted as 0; bis disulfide
  trace everywhere it appears): IMP 13 + 23 + 3 = 39; THI 0 + 19 + 4 = 23; CYS 3 + 13 + 2 = 18; blank ≈ 0
  (all three tr). So THI ≈ CYS (23 vs 18) on the pool, at ≈ 1/100 the molar dose.
- **Per mole added** (pool ng per µmol added, using the 0.0030–0.0038 mmol thiamine and 0.324 mmol cysteine
  above): thiamine 23 / 3.0–3.8 ≈ 6–8 ng/µmol; cysteine 18 / 324 ≈ 0.06 ng/µmol. Thiamine is about
  100× more effective per mole on this pool. Very rough: TIC response factor 1, technical replicates.
- **FFT:** THI 207 vs blank 142 (ratio 1.46); CYS 124 vs blank 142 (0.87); IMP 250 (1.76). SDs 21–38, so
  only IMP is clearly above blank; the author calls the 2-furylmethyl family "generally unaffected".

## 3. What it means for the model

The sweep says that at 140 °C/5 min thiamine adds almost nothing to MFT and cysteine and ribose dominate.
Madruga is 140 °C but **30 min**, in sealed bottles, in beef with ~10× precursor additions:

- **For free MFT the paper agrees with the sweep:** thiamine leaves MFT at trace (= blank), cysteine
  raises it to 3 ng/100 g, and the ribose source (IMP) to 13 ng/100 g. Ranking IMP > CYS > THI ≈ BLK.
- **For the 2-methyl-3-furyl disulfides it does not:** thiamine raises methyl 2-methyl-3-furyl disulfide
  (19 vs 13 for cysteine, from a blank at trace) and the methyl trisulfide (4 vs 2), and is the only arm
  besides IMP with bis(2-methyl-3-furyl) disulfide at trace. On the summed 2-methyl-3-furyl pool thiamine
  ≈ cysteine, at about one hundredth of cysteine's molar dose. If the model's "MFT" stands for the thiol
  plus its oxidation/methanethiol products, this paper contradicts "almost nothing at 140 °C"; if it
  stands for free MFT only, it agrees.
- Caveats on the comparison itself: 30 min vs the sweep's 5 min (more time for thiamine breakdown), MFT
  partitioning into disulfides during a 2 h, 60 °C headspace sweep is uncontrolled, and the
  thiamine-vs-cysteine difference on any single compound is within the replicate noise.

**Strength of evidence: weak.** Single internal standard on TIC with response factor 1 (the author calls the
values approximate), triplicate headspace collections of what reads as one cooked portion per arm, no
statistics, one temperature and time, no pH for three of four arms, a very large IMP dose (2.7 g/100 g) whose "ten
times" basis is not printed, numbers at the 2–25 ng/100 g level a few times above the 0.2 ng/100 g trace level. Good for
the ranking signs at 140 °C, not for rates or yields.

## What it does not give

- No temperature or time series; nothing at 100 °C.
- No thiamine + cysteine or IMP + cysteine combination arms.
- No native precursor concentrations (the "ten times" basis is not printed), no salt forms.
- No measured pH except the IMP arm adjusted to 5.6.
- No independent cooking replicates evident; no significance tests.
- No calibrated (compound-specific) quantification; no odour activity values or GC-O intensities.
- No values for bis(2-furylmethyl) disulfide in the IMP and CYS arms (blank cells).
