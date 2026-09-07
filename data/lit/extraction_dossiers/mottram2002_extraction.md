# Mottram & Nobrega 2002 — EXTRACTION (cysteine + ribose / ribose 5-phosphate / IMP, buffered and unbuffered, 140 C / 30 min)
### Do the three meat forms of ribose behave alike? No: IMP is 10-100x less reactive; ribose-5-P makes norfuraneol without cysteine; ribose needs buffer.

**Source on disk:** `data/articles/mottram2002.pdf` (owner's download, 2026-09-07). Read-only extraction from
`pdftotext -layout`; Table 1 (two pages, 12 systems x ~80 rows) and Table 2 have clean text layers and are
re-typed in full. Repo status before this dossier: cited in `docs/slr_benchmark_evaluation.md`,
`docs/reference/SCIENTIFIC_REFERENCE.md` and the two matrix-benchmark protocols; NOT in
FIT_HOLDOUT_DECLARATION and no prior dossier. ⚠ The protocol docs attribute a "pH 5, 95 C, 4 h, [13C5]ribose"
experiment and an "induction period" to this paper; neither is in it (those are Cerny & Davidek 2003/2004).
This paper is 140 C / 30 min, unlabelled, no time series.

## 0. Identity

| field | value |
|---|---|
| Title | "Formation of Sulfur Aroma Compounds in Reaction Mixtures Containing Cysteine and Three Different Forms of Ribose" |
| Authors | Donald S. Mottram and Ian C. C. Nobrega (Reading) |
| Venue | J. Agric. Food Chem. 2002, 50, 4080-4086 |
| DOI | 10.1021/jf0200826 |

## 1. Why it matters

The model's pentose lane is charged as a generic aldopentose. In meat the pentose is mostly bound as
inosine 5'-monophosphate (IMP), with smaller amounts of ribose 5-phosphate (R5P) and free ribose. This
paper heats cysteine with each of the three at 1:1, 50 mM each, under four buffer conditions, and reads
MFT, FFT, all three mercaptoketones, furfural and ~75 other volatiles. It gives (i) the reactivity
ranking IMP << R5P ~ ribose, (ii) a pH 4.2 vs 5.6 pair for every thiol in the same buffer, (iii) a
buffered-vs-unbuffered pair showing that phosphate/phthalate catalysis matters more than pH for the
2,3-enolisation (norfuraneol, MFT, mercaptoketone) branch and not for the 1,2-enolisation (furfural,
FFT) branch, and (iv) an absolute-ish norfuraneol yield from R5P alone (Table 2). It has no time axis
and no thiol-removal information.

## 2. Methods as they matter to a model

- **Charges:** separate 0.1 M solutions of cysteine and of the ribose compound (IMP disodium salt, D-ribose
  5-phosphate disodium salt, or D-ribose) in water or buffer, pH adjusted to 5.6 or 4.2 with NaOH/HCl,
  mixed in equal volumes; **6.0 mL of the mixture (0.3 mmol of each reactant) per ampule** -> **50 mmol/L
  cysteine, 50 mmol/L sugar compound, 1:1**.
- **Buffers:** phosphate = **0.2 M PYROPHOSPHATE** (disodium + tetrasodium pyrophosphate mixed to pH); phthalate
  = 0.3 M potassium hydrogen phthalate + NaOH to pH 5.6 ("0.3 M because this gave better pH control");
  unbuffered = glass-distilled water at initial pH 5.6. Conditions: unbuffered pH 5.6; phosphate pH 5.6;
  phthalate pH 5.6; phosphate pH 4.2. ⚠ Not orthophosphate — pyrophosphate; catalysis and ionic strength
  are not those of the 0.5 M KH2PO4/K2HPO4 pots the model is fitted on.
- **Vessel / heating:** 10 mL round-bottom thick-wall Pyrex ampules, flame-sealed; **CERTOclav autoclave,
  140 C, 30 min, 0.28 MPa**. Atmosphere: sealed air headspace above 6 mL.
- **pH measured before and after heating** (rows in Table 1). Unbuffered: R5P 5.6 -> 4.2, ribose 5.6 -> 3.9,
  IMP 5.6 -> 5.7. Buffers held within 0.3 unit (phosphate 5.6 -> 5.5; phthalate 5.6 -> 5.3-5.4; phosphate
  4.2 -> 3.9, IMP 4.5).
- **Replicates: triplicate** per system; Table 1 = means. "The mean coefficient of variance (CV) for
  quantities of individual components was 22% and, with the exception of some compounds that were present in
  relatively small amounts, no compound showed a CV >40%."
- **Isolation:** reaction mixture + 20 mL of the same buffer (or water) in a 250 mL flask at 60 C, swept with
  N2 at 40 mL/min for 1.5 h onto 85 mg Tenax-TA; internal standard 1,2-dichlorobenzene 130 ng added to the
  trap afterwards; purged 5 min.
- **GC-MS:** HP 5890/5972, BPX5 50 m x 0.32 mm, thermal desorption 250 C, cryofocus 0 C, 4 C/min to 280 C;
  EI 70 eV. Some samples also on BP20 for LRI.
- **Quantification: RELATIVE.** "Approximate quantities of the volatiles in the concentrated headspace ...
  estimated by comparison of their peak areas, in the total ion current chromatogram, with that of the
  1,2-dichlorobenzene internal standard using a response factor of 1. This allowed comparison of the relative
  contributions the volatiles made to the headspaces of the different systems but did not provide absolute
  concentrations in the aqueous solutions." Unit: **ng per mmol of sugar** (footnote a) — headspace-recovered
  amount per 0.3 mmol sugar charged. tr = < 0.5 ng/mmol; "-" = below detection (~0.1 ng/mmol); "+" =
  present but confounded by an adjacent peak.
- **Sugar-only blanks (Table 2):** 0.3 mmol R5P or ribose in 6 mL water or phosphate buffer, heated the same
  way, extracted with 6 mL dichloromethane, dichlorobenzene IS (36 µg) added, concentrated, 1 µL splitless;
  "quantified using the internal standard" — µg per mmol sugar, mean +- SD of triplicates. Still IS-equivalents
  (no response-factor statement), but a solvent extraction of the whole solution, so closer to a solution
  yield than the headspace numbers.

## 3. Tables re-typed

### Table 1. "Approximate Quantities(a) of Volatiles Identified in the Headspace of Heated Cysteine Model Systems Containing Ribose 5-Phosphate (Rib-PO4), Ribose, or Inosine 5'-Monophosphate (IMP)(b)"

Footnote a: "Approximate quantities in headspace (ng/mmol of sugar) given as means of triplicate analyses; tr,
trace (<0.5 ng/mmol of sugar); -, below detection limit (~0.1 ng/mmol of sugar); +, present in small amounts
and quantification confounded by adjacent peak." b: "Each model system consisted of 0.3 mmol of cysteine and
0.3 mmol of the ribose-containing compound in 6 mL of water or buffer." c: MS and LRI agree with authentic
samples. d: pair of diastereoisomers. MS-ref column omitted.

Column order in every block: **U-R5P, U-rib, U-IMP | P5.6-R5P, P5.6-rib, P5.6-IMP | Ph5.6-R5P, Ph5.6-rib,
Ph5.6-IMP | P4.2-R5P, P4.2-rib, P4.2-IMP** (U = unbuffered, P = phosphate, Ph = phthalate).

| row | U-R5P | U-rib | U-IMP | P5.6-R5P | P5.6-rib | P5.6-IMP | Ph5.6-R5P | Ph5.6-rib | Ph5.6-IMP | P4.2-R5P | P4.2-rib | P4.2-IMP |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| initial pH | 5.6 | 5.6 | 5.6 | 5.6 | 5.6 | 5.6 | 5.6 | 5.6 | 5.6 | 4.2 | 4.2 | 4.2 |
| final pH | 4.2 | 3.9 | 5.7 | 5.5 | 5.5 | 5.7 | 5.4 | 5.3 | 5.4 | 3.9 | 3.9 | 4.5 |
| **3-mercapto-2-butanone** | 2670 | 29 | 10 | 8860 | 9820 | 43 | 5570 | 3520 | 70 | 3890 | 2090 | 500 |
| **3-mercapto-2-pentanone** | 9220 | 850 | tr | 2490 | 4020 | 23 | 2370 | 1540 | 40 | 9260 | 6600 | 700 |
| **2-mercapto-3-pentanone** | 1870 | 40 | 7 | 2620 | 3280 | 37 | 1640 | 921 | 53 | 2200 | 1330 | 533 |
| total mercaptoketones | 13760 | 919 | 17 | 13970 | 17120 | 103 | 9580 | 5981 | 163 | 15350 | 10020 | 1733 |
| **2-methyl-3-furanthiol** | 2940 | 1010 | 97 | 2050 | 2230 | 67 | 2960 | 2620 | 143 | 10230 | 8250 | 1867 |
| **2-furanmethanethiol** (FFT) | 3310 | 2610 | 13 | 496 | 1750 | 23 | 840 | 3540 | 47 | 2360 | 3160 | 5000 |
| 3-thiophenethiol | 2400 | 37 | 50 | 716 | 476 | 123 | 1100 | 353 | 213 | 3770 | 848 | 367 |
| 2-methyl-3-thiophenethiol | 699 | 55 | 20 | 1260 | 730 | 23 | 1180 | 663 | 40 | 1390 | 1050 | 400 |
| 2-thiophenemethanethiol | 104 | 16 | - | 124 | 36 | tr | 171 | 67 | tr | 128 | 66 | 21 |
| total furan- and thiophenethiols | 9453 | 3728 | 180 | 4646 | 5222 | 236 | 6251 | 7243 | 443 | 17878 | 13374 | 7655 |
| bis(1-methyl-2-oxopropyl) disulfide (d) | 20 | - | - | 262 | 131 | - | 458 | 62 | - | tr | tr | tr |
| 3-(1-methyl-2-oxopropyldithio)pentan-2-one (d) | 212 | tr | - | 181 | 10 | - | 190 | 95 | - | 9 | tr | tr |
| 2-(1-methyl-2-oxopropyldithio)pentan-3-one (d) | 16 | tr | - | 202 | 200 | - | 233 | 66 | - | tr | tr | tr |
| bis(1-ethyl-2-oxopropyl) disulfide (d) | 333 | 34 | - | 42 | 29 | - | 49 | 23 | - | 24 | 12 | tr |
| 3-(1-methyl-2-oxobutyldithio)pentan-2-one (d) | 111 | 7 | - | 60 | 53 | - | 47 | 32 | - | 4 | 9 | tr |
| bis(1-methyl-2-oxobutyl) disulfide (d) | - | - | - | 40 | 33 | - | 12 | tr | - | - | - | tr |
| total oxoalkyl disulfides | 692 | 41 | - | 787 | 456 | - | 989 | 278 | - | 37 | 21 | tr |
| 1-[2-methyl-(3-furyldithio)]propan-2-one | 41 | 13 | - | - | - | - | - | - | - | 12 | 27 | 3 |
| 2-[2-methyl-(3-furyldithio)]butan-3-one | 71 | 5 | tr | 146 | 149 | tr | 580 | - | - | 43 | 42 | 17 |
| **bis(2-methyl-3-furyl) disulfide** | 77 | 56 | 7 | 33 | 87 | tr | 240 | 55 | tr | 167 | 319 | 63 |
| 3-[(2-methyl-(3-furyldithio)]pentan-2-one | 281 | 38 | - | 86 | 35 | - | 180 | 117 | - | 157 | 201 | 30 |
| 2-[(2-methyl-(3-furyldithio)]pentan-3-one | 48 | 8 | tr | + | 74 | - | 152 | 83 | - | 34 | 37 | 13 |
| 2-(2-furylmethyldithio)butan-3-one | 33 | 3 | - | 18 | 112 | - | 53 | 353 | - | 3 | 3 | tr |
| 2-methyl-3-(2-furylmethyldithio)furan | 47 | 47 | - | + | + | - | 18 | 148 | - | 14 | 32 | tr |
| 3-(2-furylmethyldithio)pentan-2-one | 215 | 63 | - | tr | 35 | - | 37 | 102 | - | 11 | 22 | tr |
| 2-(2-furylmethyldithio)pentan-3-one | + | + | - | 13 | 52 | - | 12 | 87 | - | 3 | 3 | tr |
| bis(2-furylmethyl) disulfide | 43 | 80 | - | - | 20 | - | - | 68 | - | tr | 7 | tr |
| total furyl disulfides | 856 | 313 | 7 | 296 | 564 | - | 1272 | 1013 | - | 444 | 693 | 126 |
| 2-(3-thienyldithio)butan-3-one | + | - | - | 94 | 11 | tr | 61 | tr | tr | 11 | tr | - |
| 2-methyl-3-(3-thienyldithio)furan | 80 | tr | 7 | 38 | 14 | tr | 54 | 45 | 7 | 63 | 17 | 13 |
| 2-methyl-3-[2-methyl-(3-thienyldithio)]furan | 49 | 11 | tr | 63 | 23 | tr | 88 | 45 | tr | 24 | 31 | 17 |
| 2-(3-thienyldithio)pentan-3-one | 46 | tr | - | 27 | 8 | 10 | 14 | tr | - | 8 | tr | tr |
| 3-(2-furylmethyldithio)thiophene | + | 7 | - | - | - | - | - | - | - | - | - | - |
| bis(3-thienyl) disulfide | 49 | tr | 3 | - | - | 20 | - | - | tr | 8 | tr | tr |
| 2-methyl-3-(3-thienyldithio)thiophene | 25 | tr | tr | 13 | tr | 10 | tr | tr | tr | 8 | tr | tr |
| bis(2-methyl-3-thienyl) disulfide | 9 | tr | tr | 21 | tr | - | tr | tr | - | tr | tr | 3 |
| total thienyl disulfides | 258 | 18 | 10 | 256 | 56 | 40 | 217 | 90 | 7 | 122 | 48 | 33 |
| 2-methylthiophene | 1220 | 300 | 70 | 744 | 550 | 33 | 1180 | 537 | 73 | 1470 | 828 | 833 |
| 4,5-dihydro-2-methylthiophene | - | - | - | 123 | 199 | - | 403 | 246 | - | - | - | + |
| 2-ethylthiophene | - | - | - | 51 | 144 | 3 | 180 | 50 | 3 | - | - | - |
| 2,3-dimethylthiophene | 62 | tr | tr | 96 | 201 | tr | 261 | 71 | tr | - | - | - |
| 2-formylthiophene | 146 | 68 | tr | 133 | 258 | tr | 239 | 239 | tr | 238 | 310 | 20 |
| 5-methyl-2-formylthiophene | 55 | 3 | 3 | 27 | 17 | 13 | 78 | 18 | 3 | tr | tr | tr |
| 3-methyl-2-formylthiophene | 190 | 8 | 7 | 226 | 359 | 10 | 500 | 462 | 7 | 110 | 173 | 17 |
| 2-acetyl-3-methylthiophene | 46 | - | - | 54 | - | - | 60 | tr | - | 79 | 6 | tr |
| 2-propanoylthiophene | 17 | - | - | 110 | tr | tr | 152 | tr | tr | 24 | 32 | tr |
| 3-ethyl-2-formylthiophene | 60 | 8 | tr | 69 | 112 | tr | 118 | 132 | tr | 111 | 266 | 20 |
| dimethylformylthiophene | 291 | 20 | - | 63 | 26 | tr | 137 | 44 | tr | 132 | 129 | 20 |
| total thiophenes | 2087 | 407 | 80 | 1696 | 1866 | 59 | 3308 | 1799 | 86 | 2164 | 1744 | 910 |
| (E or Z)-3,5-dimethyl-1,2-dithiolan-4-one | 457 | 9 | 27 | 224 | 116 | 20 | 361 | 119 | 30 | 1224 | 113 | 297 |
| (E or Z)-3,5-dimethyl-1,2-dithiolan-4-one | 490 | 7 | 20 | 131 | 72 | 13 | 240 | 70 | 23 | 1019 | 209 | 237 |
| 3-ethyl-1,2-dithiolan-4-one | 21 | - | - | 36 | 22 | - | 44 | 20 | - | 24 | 24 | 7 |
| 3-methyl-1,2-dithian-4-one | 190 | 6 | 3 | 128 | 307 | 10 | 222 | 453 | 10 | 290 | 458 | 47 |
| total dithianones and dithiolanones | 1158 | 22 | 50 | 519 | 517 | 43 | 867 | 662 | 63 | 2557 | 804 | 588 |
| 4,5-dihydro-3(2H)-thiophenone | 28 | tr | 3 | 76 | 147 | 17 | 108 | 90 | 10 | 33 | 23 | 7 |
| 4,5-dihydro-5-methyl-3(2H)-thiophenone | 117 | 17 | tr | + | + | + | + | + | + | 251 | 51 | 20 |
| 4,5-dihydro-2-methyl-3(2H)-thiophenone | 244 | 26 | 10 | 3040 | 1890 | 23 | 2200 | 1020 | 33 | 868 | 536 | 87 |
| dihydro-2,(4 or 5)-dimethyl-3(2H)-thiophenone | tr | - | - | 189 | 159 | tr | 253 | 104 | tr | 17 | 22 | tr |
| dihydro-2,(4 or 5)-dimethyl-3(2H)-thiophenone | tr | - | - | 306 | 260 | 10 | 302 | 123 | tr | 22 | 26 | tr |
| dihydro-2 or 5-ethyl-3(2H)-thiophenone | - | - | - | 138 | 90 | tr | 201 | 81 | tr | - | - | tr |
| ethyl-3(2H)-thiophenone | - | - | - | 40 | 15 | - | 52 | 12 | - | - | - | - |
| total thiophenones | 389 | 43 | 13 | 3789 | 2561 | 50 | 3116 | 1430 | 43 | 1191 | 658 | 114 |
| 2,3-dihydro-6-methylthiothieno[2,3-c]furan | tr | - | - | 224 | 1320 | - | 406 | 840 | - | tr | 43 | - |
| thieno[2,3-b]thiophene | 290 | 8 | 10 | + | 43 | 23 | 147 | 47 | 23 | 275 | 37 | 3 |
| thieno[3,2-b]thiophene | 7 | tr | tr | 23 | 24 | 7 | 39 | 60 | tr | - | - | - |
| a dihydrothienothiophene | 68 | tr | - | 887 | 1710 | 7 | 1310 | 1490 | 3 | 137 | 96 | 7 |
| a methyldihydrothienothiophene | 5 | tr | - | 139 | 628 | - | 243 | 660 | - | tr | tr | tr |
| a methyldihydrothienothiophene | 73 | 3 | - | 156 | 110 | - | 400 | 143 | - | 33 | 26 | tr |
| a methyldihydrothienothiophene | 153 | 5 | - | 140 | 103 | - | 334 | 128 | - | 68 | 58 | tr |
| a dimethyldihydrothienothiophene | 30 | tr | - | 152 | 254 | - | 497 | 460 | - | 7 | 16 | tr |
| total bicyclic compounds | 626 | 16 | 10 | 1721 | 4192 | 37 | 3376 | 3828 | 26 | 520 | 276 | 10 |
| (E)-3,5-dimethyl-1,2,4-trithiolane | - | - | 7 | tr | tr | 23 | + | tr | tr | - | - | - |
| (Z)-3,5-dimethyl-1,2,4-trithiolane | - | - | tr | tr | tr | 20 | + | tr | tr | - | - | - |
| 3-methyl-1,2,4-trithiane | - | - | 7 | 27 | 8 | 97 | 38 | 6 | 7 | - | - | - |
| 1,2,4,5-tetrathiane | - | - | 63 | - | - | 210 | - | - | 67 | - | - | - |
| total trithiolanes, trithianes, and tetrathianes | - | - | 77 | 27 | 8 | 350 | 38 | 6 | 74 | - | - | - |
| 2-pentanone | 124 | 12 | 27 | 801 | 388 | 40 | 1040 | 469 | 67 | 352 | 123 | 117 |
| 3-pentanone | 26 | - | 13 | 457 | 388 | 33 | 703 | 279 | 43 | 82 | 22 | 37 |
| 2,3-pentanedione | + | + | - | + | + | + | + | + | + | + | + | + |
| 3-hydroxy-2-butanone | tr | - | - | 453 | 98 | - | 1060 | tr | - | + | + | - |
| 2,4-pentanedione | 2530 | 10 | 3 | 912 | 102 | 13 | 1170 | 33 | 13 | 5640 | 482 | 367 |
| methylpyrazine | - | - | - | 81 | 267 | - | tr | 294 | - | tr | tr | - |
| **2-furfural** | **1460** | **4990** | **tr** | **-** | **tr** | **tr** | **-** | **260** | **tr** | **1260** | **2290** | **27** |
| total non-sulfur compounds | 4140 | 5012 | 43 | 2704 | 1243 | 86 | 3973 | 1335 | 123 | 7334 | 2917 | 665 |

### Table 2. "Approximate Quantities(a) (Micrograms per Millimole of Sugar) of Major Compounds Identified in Dichloromethane Extracts of Heated Ribose 5-Phosphate or Ribose Solutions(b)"

Footnote a: "Quantities are the mean of triplicate analyses with the standard deviation in parentheses; tr,
trace (<5 µg/mmol); -, not detected." b: "Each solution consisted of 0.3 mmol of ribose or ribose 5-phosphate
in 6 mL of water or phosphate buffer." NO CYSTEINE in these runs.

| compound | unbuffered Rib-PO4 | unbuffered ribose | phosphate Rib-PO4 | phosphate ribose |
|---|---:|---:|---:|---:|
| final pH (initial pH 5.6) | 4.4 | 3.8 | 5.7 | 5.7 |
| 2-furfural | 263 (63) | 103 (23) | tr | 107 (3) |
| 4-hydroxy-5-methyl-3(2H)-furanone (norfuraneol) | 280 (27) | - | 287 (47) | tr |

In mol % of sugar charged (MW furfural 96.09, norfuraneol 114.10): furfural from R5P unbuffered **0.27 %**,
from ribose unbuffered **0.107 %**, from ribose in phosphate **0.111 %**; norfuraneol from R5P **0.245 %**
(unbuffered, final pH 4.4) and **0.252 %** (phosphate, pH 5.7); from ribose < 0.004 % (tr) or nd. Text:
norfuraneol "could not be detected in the headspace of the heated reaction mixtures" (Table 1 systems) —
it is water-soluble and not Tenax-recoverable at 60 C.

## 4. What the repo could take

All Table 1 numbers are headspace ng per mmol sugar, response factor 1, CV ~22 %: use as ratios only.
Table 2 is closer to a solution yield.

### 4.1 Candidate near-absolute rows (Table 2, sugar alone, 140 C / 30 min, 50 mM, no cysteine)

| system | product | µg/mmol sugar (SD) | mol % | comment |
|---|---|---:|---:|---|
| ribose 50 mM, 0.2 M pyrophosphate pH 5.6 -> 5.7 | furfural | 107 (3) | **0.11** | unbuffered (pH -> 3.8) gives the same 103 (23): furfural yield from ribose alone is pH-flat between 3.8 and 5.7 at this severity |
| R5P 50 mM, pyrophosphate pH 5.7 | norfuraneol | 287 (47) | **0.25** | R5P dephosphorylates to NF without any amine; ribose gives tr — the sugar-only NF route is a phosphate-ester route, not a caramelisation one |
| R5P 50 mM, unbuffered -> pH 4.4 | norfuraneol / furfural | 280 (27) / 263 (63) | 0.25 / 0.27 | |
| ribose 50 mM, either | norfuraneol | tr / - | < 0.004 | ribose alone makes no NF in 30 min at 140 C |

Caveat: IS-equivalents (dichlorobenzene), not calibrated; DCM extraction efficiency for NF (very polar) is
unknown, so the NF numbers are lower bounds if anything.

### 4.2 Within-study ratios on the thiols (Table 1)

**pH 4.2 vs 5.6, same pyrophosphate buffer, 140 C / 30 min** (ratio = P4.2 / P5.6):

| compound | R5P | ribose | IMP |
|---|---:|---:|---:|
| MFT | 10230 / 2050 = **5.0** | 8250 / 2230 = **3.7** | 1867 / 67 = 28 |
| FFT | 2360 / 496 = 4.8 | 3160 / 1750 = **1.8** | 5000 / 23 = 217 |
| 3-mercapto-2-pentanone | 9260 / 2490 = 3.7 | 6600 / 4020 = **1.6** | 700 / 23 = 30 |
| 2-mercapto-3-pentanone | 2200 / 2620 = 0.84 | 1330 / 3280 = **0.41** | 533 / 37 = 14 |
| 3-mercapto-2-butanone | 3890 / 8860 = 0.44 | 2090 / 9820 = **0.21** | 500 / 43 = 12 |
| furfural | 1260 / - (> 10000) | 2290 / tr (> 4000) | 27 / tr |
| bis(2-methyl-3-furyl) disulfide | 167 / 33 = 5.1 | 319 / 87 = 3.7 | 63 / tr |
| dithiolanones + dithianones | 2557 / 519 = 4.9 | 804 / 517 = 1.6 | 588 / 43 = 14 |
| thiophenones | 1191 / 3789 = 0.31 | 658 / 2561 = 0.26 | 114 / 50 = 2.3 |
| MFT / FFT | 4.3 (pH 4.2) vs 4.1 (5.6) | **2.6 (4.2) vs 1.27 (5.6)** | 0.37 vs 2.9 |

Directional: **MFT and furfural rise as pH falls from 5.6 to 4.2; 2-mercapto-3-pentanone, 3-mercapto-2-
butanone and the thiophenones fall.** The two mercaptopentanone isomers again move in opposite directions
with pH (cf. Whitfield 2001, Cerny 2007). Furfural goes from trace at pH 5.5 to 2290 at pH 3.9 — the same
switch Cerny 2007 sees between pH 5.5 and 6 at 145 C, seen here from the acid side.

**Ribose vs ribose 5-phosphate vs IMP, same buffer (ratio to ribose):**

| condition | MFT R5P : rib : IMP | FFT R5P : rib : IMP | 3-MP R5P : rib : IMP | furfural R5P : rib : IMP |
|---|---|---|---|---|
| unbuffered (final pH 4.2 / 3.9 / 5.7) | 2.9 : 1 : 0.10 | 1.27 : 1 : 0.005 | 10.8 : 1 : ~0 | 0.29 : 1 : ~0 |
| phosphate 5.6 | 0.92 : 1 : 0.03 | 0.28 : 1 : 0.013 | 0.62 : 1 : 0.006 | - : tr : tr |
| phthalate 5.6 | 1.13 : 1 : 0.055 | 0.24 : 1 : 0.013 | 1.54 : 1 : 0.026 | - : 260 : tr |
| phosphate 4.2 | 1.24 : 1 : 0.23 | 0.75 : 1 : 1.58 | 1.40 : 1 : 0.11 | 0.55 : 1 : 0.012 |

Directional: **buffered at pH 5.6, ribose and R5P give the same MFT within 10-25 % and the same 3-MP within
~1.6x; R5P gives 3.5-4x LESS FFT and no furfural** (R5P bypasses the Amadori / 3-deoxyosone route that makes
furfural). **IMP is 10-100x less reactive at pH 5.6** ("most compounds were found in the headspace at
concentrations 10-100 times lower"), rising to ~0.1-0.25x of ribose at pH 4.2 (glycoside hydrolysis is
acid-catalysed) — with the striking exception that **IMP at pH 4.2 gives the MOST FFT of any system (5000)**
and a MFT/FFT of 0.37. Unbuffered, ribose is "relatively unreactive" (MFT 1010, 3-MP 850) while R5P is not
(2940, 9220).

**Buffer effect on ribose + cysteine (ratio to unbuffered; note unbuffered ended at pH 3.9):**

| compound | phosphate 5.6 / U | phthalate 5.6 / U | phosphate 4.2 / U |
|---|---:|---:|---:|
| MFT | 2230 / 1010 = 2.2 | 2620 / 1010 = 2.6 | 8250 / 1010 = **8.2** |
| 3-mercapto-2-pentanone | 4.7 | 1.8 | 7.8 |
| 2-mercapto-3-pentanone | 82 | 23 | 33 |
| 3-mercapto-2-butanone | 339 | 121 | 72 |
| FFT | 0.67 | 1.36 | 1.21 |
| furfural | tr / 4990 (< 0.0001) | 0.052 | 0.46 |
| thiophenones | 60 | 33 | 15 |

Authors: "the marked increase in these compounds in the buffered ribose systems, at both pH 5.6 and 4.2,
compared with the unbuffered ribose system, suggests that catalysis by the buffer occurred, which was
greater than any effect of pH"; "2-furylmethanethiol, which is formed from 2-furfural, also showed little or
no increase in the presence of buffer"; "both the phosphate and phthalate buffers exhibited similar extents
of catalysis, indicating that acid-base catalysis was the dominant mechanism." For the model: **the buffer
catalyses the 2,3-enolisation branch (NF / 1-deoxyosone -> MFT, mercaptoketones) by 2-8x and the
1,2-enolisation branch (furfural -> FFT) not at all.** The comparison is confounded by the unbuffered pH
drop to 3.9, which by the pH-4.2 column should have RAISED MFT; that the buffered pH-5.6 pot still beats
it 2.2x is the catalysis signal.

### 4.3 Other numbers worth a line

- The unbuffered systems show the acid drift the model's pH-trajectory wave cares about: ribose + cysteine
  5.6 -> 3.9 and R5P + cysteine 5.6 -> 4.2 in 30 min at 140 C; IMP + cysteine 5.6 -> 5.7 (no reaction, no
  acid). Sugar alone (Table 2): R5P 5.6 -> 4.4 (phosphate release), ribose 5.6 -> 3.8.
- Disulfides are 3-10 % of their parent thiols everywhere (e.g. bis(2-methyl-3-furyl) disulfide / MFT =
  0.016-0.04 in phosphate; oxoalkyl disulfides / mercaptoketones = 0.03-0.10): under a 30 min sealed
  autoclave, oxidative dimerisation is a minor sink.
- Trithiolanes / trithianes / tetrathiane (cysteine self-degradation products) appear only where the sugar
  does not react (IMP at pH 5.6: 350) — the authors read this as competition for cysteine-breakdown
  intermediates.

### 4.4 Suggested role

Not previously declared. Everything is relative and at 140 C / 30 min in pyrophosphate; suitable as
**directional / ratio HOLD-OUT** (the pH 4.2 vs 5.6 ratios in §4.2 and the R5P ~ ribose equivalence),
not as levels. The Table 2 NF and furfural yields from sugar alone are the only near-absolute numbers and
could serve as VALIDATE rows for a sugar-only run if the model is ever charged with R5P.

## 5. Caveats

1. **Relative headspace quantities, RF = 1, dynamic headspace at 60 C on Tenax.** Not solution
   concentrations; the authors say so explicitly.
2. **Pyrophosphate buffer (0.2 M), not orthophosphate.** Buffer catalysis is shown to matter 2-8x, so these
   levels are not commensurable with 0.5 M KH2PO4 pots even in ratio form across studies.
3. **Unbuffered runs drift 1.7 pH units**; "unbuffered pH 5.6" is really a pH 5.6 -> 3.9 trajectory.
4. **Single time point (30 min)**; no thiol-removal information beyond the disulfide/thiol ratios.
5. IMP was charged as the disodium salt and R5P as the disodium salt; ionic strength differs from the ribose
   systems even unbuffered.
6. Norfuraneol not detectable by the headspace method; its Table 2 numbers come from a different workup
   (DCM, no cysteine) and cannot be paired with the Table 1 thiols.
7. The repo's protocol docs mis-cite this paper for a 95 C / 4 h / 13C5 experiment and an "induction
   period"; those belong to Cerny & Davidek 2003/2004.
