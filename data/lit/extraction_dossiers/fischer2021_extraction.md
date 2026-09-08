# Fischer, Cachon & Cayot 2021 — EXTRACTION (HS-SPME GC-MS of a commercial pea protein isolate at extraction pH 6.5 / 4.5 / 2.0; nine beany volatiles semi-quantified in ug per g isolate)
### The off-flavour a pea isolate carries into a recipe, measured as free-plus-releasable volatile per gram of powder — a binding/release study, not a lipoxygenase study.

**Source on disk:** `data/articles/Fischer2021.pdf` (owner's download, 2026-09-08; Food Research
International 150 (2021) 110760, 9 pages). Read from the scratchpad text layer
(`Fischer2021.txt`); Tables 1, 2 and 4 are clean; Table 3 is clean except one cell (2-ethylfuran at
pH 2, garbled in the PDF's own text layer — re-extracted with pypdf in plain and layout mode,
still unreadable, marked UNREADABLE below). Fig. 1 (siloxane vs number of extractions) and Fig. 2
(binding diagram) are FIGURE-ONLY and carry nothing the repository needs. Repo status before this
dossier: LOX-01 holds raw pea MILK levels (Zhang 2020b) and one isolate LOX residue (Gao 2020) but
no volatile level for a finished pea isolate; Gao 2020's isolate hexanal is FIGURE-ONLY.

## 0. Identity

| field | value |
|---|---|
| Title | "Effects of extraction pH on the volatile compounds from pea protein isolate: Semi-Quantification method using HS-SPME-GC-MS" |
| Authors | Estelle Fischer, Remy Cachon, Nathalie Cayot (corresponding) — Univ. Bourgogne Franche-Comte, AgroSup Dijon, PAM UMR A 02.102 |
| Venue | Food Research International 150 (2021) 110760; received 8 Feb 2021, accepted 13 Oct 2021, online 17 Oct 2021 |
| DOI | 10.1016/j.foodres.2021.110760 |
| Funding | Regional Council Bourgogne-Franche-Comte, ERDF, and a grant from Roquette Freres S.A. (the isolate's manufacturer) |
| Companion | Fischer, Cachon & Cayot 2020, Trends Food Sci Technol 95, 196 (pea vs soy review, cited by Bi 2026 for the citric-acid quench) |
| "Extraction pH" | means the pH of the 2 mL suspension in the SPME vial, i.e. the pH at which volatiles are RELEASED from the powder for headspace sampling — not the pH at which the isolate was manufactured |

## 1. Why it matters

The repository needs the beany load that a pea protein isolate brings with it before any heating
(the "carried-over off-flavour"), in a unit per gram of isolate that a recipe can be charged with.
This paper gives exactly that for hexanal, nonanal, (E)-2-nonenal, 3-methylbutanal, benzaldehyde,
1-octen-3-ol, 3-octen-2-one, 2-pentylfuran and 2,5-dimethylpyrazine in two lots of one commercial
spray-dried isolate (Roquette, 85 % protein), with an explicit external calibration, LOD/LOQ and
repeatability. It also shows how much of each compound is protein-bound at neutral pH and released
by acidification (hexanal +59 % at pH 4.5), which is the binding side of the same story the repo's
`protein_matrices.yml` tries to carry. It contains **no lipoxygenase measurement and no time
course**; its numbers are levels and within-study pH ratios.

## 2. Methods as they matter to a model

- **Isolate:** "Spray-dried pea protein isolate (85 % protein dry matter, composed mainly of
  globulins) was supplied by Roquette Freres S.A." Two samples, PPI1 and PPI2, of unstated
  provenance / age ("It could be hypothesized that the two samples of PPI do not have a similar
  age"). No lipid content, moisture, LOX activity or manufacturing pH is given.
- **Suspension:** 0.2 g PPI weighed into a clear 20 mL vial; liquid added to **2 mL, 10 % (w/v)**;
  liquid/gas 2/18 (v/v). "Distilled water was used for extraction at neutral pH (6.5), 0.1 M HCl
  for extraction at pH 4.5, and 1 M HCl ... for extraction at pH 2.0." Triplicate per pH. The pH
  values are the resulting suspension pHs as the authors label them; no buffer.
- **Rationale for the three pHs:** 6.5 = "neutral, non-denaturing pH, corresponding to the pH of
  the protein matrix"; 4.5 = "partially denaturing (7S globulin precipitation)"; 2.0 = "strong
  denaturing ... (including 11S globulins)".
- **HS-SPME:** DVB/CAR/PDMS 50/30 um, 1 cm, manual holder; **40 C**, 350 rpm, dark; **equilibration
  30 min + extraction 60 min** (90 min total at 40 C in the acid or water). Fibre conditioned 270 C
  / 30 min; desorption 250 C / 5 min, split mode with 140 mL/min purge at 0 min.
- **GC-MS:** HP 6890 / HP 5973 quadrupole; DB-WAX 30 m x 0.32 mm x 0.25 um; He 1.4 mL/min (43
  cm/s); 40 C (3 min) -> 100 C at 3 C/min -> 230 C at 5 C/min (10 min); 59 min run; source 230 C,
  transfer line 190 C; EI 70 eV, full scan m/z 29-400; NIST 08 / Wiley / INRA libraries;
  **integration limit 50 000 area counts**; the nine target compounds confirmed against standards.
- **Calibration solution:** ~50 mg (exactly weighed) of each of the nine compounds in 1 L distilled
  water (50 ppm each), overnight homogenisation; 20 uL into 2 mL water in the vial = **0.5 ppm
  (0.5 mg/L)** each for the repeatability test (triplicate, same day, two fibres one month apart).
- **Semi-quantification, verbatim:** "Naphthalene D8 ... was first used as an internal standard,
  but it gave poor repeatability. ... an external calibration was chosen. The calibration curves of
  each nine compounds of interest were obtained for concentrations ranging from 0.001 to 2.5 ppm,
  in distilled water. Each compound was analyzed in the presence of the others in order to take
  into account the potential competition between compounds in the headspace. For two compounds,
  hexanal and 3-methylbutanal, the calibration curves were no longer linear above a given
  concentration. Consequently, the two compounds were run separately to obtain the calibration
  curves." And: "this method did not take into account the extraction yield, as internal
  calibration would have done, and was thus called a semi-quantification method." So: **external
  calibration in pure water, no internal standard, no matrix matching, no correction for the
  isolate's own headspace partition (protein binding lowers recovery); result expressed as ug of
  compound per g of isolate powder.**
- **Worked example, verbatim (PPI1, pH 2.0, hexanal):** area 7 078 843 A.U.; [hexanal]assay (ug/mL)
  = (7 078 843 - 157 037) / 1 x 10^7 = 0.69 ug/mL; [hexanal]sample (ug/g) = 0.69 x 2 mL / 0.2041 g
  = 6.78 ug/g; mean of n = 3: 6.7 +/- 0.2 ug hexanal / g isolate. (Matches Table 4 PPI1 pH 2.0.)
- **Table 3 belongs to PPI1** (not stated in the caption; verified here: Table 3 area/g x 0.2 g,
  through the Table 2 hexanal, nonanal, benzaldehyde and 2-pentylfuran equations, reproduces the
  PPI1 column of Table 4 to the printed precision — 5.08, 8.07, 6.67 ug/g for hexanal).
- **Statistics:** one-way ANOVA + Tukey, p < 0.05, letters in Tables 3 and 4.
- **Method repeatability:** "variation coefficient of the method for a single fiber was 15 %"
  (abstract); 0.73-17.67 % by compound and fibre (Table 1); between fibres 6.2-67.9 %; fibre life
  ~100 extractions monitored by siloxane bleed (Fig. 1).
- **LOX statement:** none. The word lipoxygenase does not occur; the formation pathways are cited
  only to justify choosing chemically diverse calibration compounds (Frankel 1983; Gargouri 2008;
  Sessa 1979; Ullrich & Grosch 1987). ⚠ Nothing in the protocol inactivates residual LOX (Gao
  2020: a finished pea isolate keeps 18.6 U/g), and the pH 6.5 vial spends 90 min at 40 C in
  water with 10 % powder; whether any hexanal is formed during the SPME step itself is not
  controlled (flag 3).

## 3. Tables re-typed

### Table 1. "Information and chromatographic data on the studied volatile compounds. Variation coefficient for the calibration method of the two SPME fibers, at 0.5 ppm of each compound (n = 3)."

Solubility (g/L) from ChemSpider / PubChem; log P and Henry's constant (Pa m3/mol) calculated by
EPI Suite 4.1. Areas are mean of 3 at 0.5 ppm in water.

| family | compound | CAS | solubility g/L | log P | Henry Pa m3/mol | RT min | fibre 1 mean area | CV % | fibre 2 mean area | CV % | CV between fibres % |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| aldehyde | hexanal | 66-25-1 | 3.520 | 1.80 | 2.13E+01 | 5.3 | 4 373 173 | 6.18 | 4 207 065 | 4.30 | 18.85 |
| | nonanal | 124-19-6 | 0.096 | 3.27 | 5.00E+01 | 17.0 | 9 002 917 | 2.84 | 3 846 293 | 3.71 | 40.5 |
| | 2-nonenal | 18829-56-6 | 0.205 | 3.06 | 2.34E+01 | 22.5 | 4 774 558 | 12.80 | 2 866 753 | 8.79 | 28.25 |
| | 3-methylbutanal | 590-86-3 | 11.230 | 1.23 | 1.61E+01 | 2.4 | 884 178 | 14.80 | 806 879 | 3.53 | 6.18 |
| | benzaldehyde | 100-52-7 | 6.100 | 1.71 | 1.36E+00 | 21.8 | 4 262 010 | 14.71 | 3 512 612 | 1.28 | 11.14 |
| alcohol | 1-octen-3-ol | 3391-86-4 | 1.840 | 2.60 | 2.34E+00 | 19.8 | 9 121 866 | 10.45 | 7 425 840 | 1.11 | 13.15 |
| ketone | 3-octen-2-one | 18402-82-9 | 1.050 | 2.29 | 9.71E+00 | 17.5 | 8 515 285 | 14.01 | 6 493 417 | 0.73 | 16.79 |
| furan | 2-pentylfuran | 3777-69-3 | 0.041 | 3.87 | 1.87E+03 | 10.4 | 2 340 873 | 7.64 | 108 573 | 17.67 | 67.9 |
| pyrazine | 2,5-dimethylpyrazine | 123-32-0 | 31.970 | 1.03 | 3.60E-01 | 13.7 | 712 528 | 10.44 | 652 742 | 1.50 | 6.19 |

### Table 2. "Calibration curves, LOD, LOQ, linearity of the compounds studied."

y = peak area (A.U.), x = ppm (= mg/L = ug/mL) in the 2 mL assay, pure water, 40 C, 30 + 60 min.
LOD 0 = detected with no area limit; LOD 50 k = detected above the 50 000 integration limit; LOQ =
lowest concentration quantifiable within the linear range. Slopes and intercepts as printed
(one significant figure on the slopes).

| compound | LOD 0 ppm | LOD 50 k ppm | LOQ ppm | linearity ppm | equation | R2 |
|---|---:|---:|---:|---|---|---:|
| hexanal | < 0.001 | < 0.001 | 0.05 | 0.05-2.5 | y = 1 x 10^7 x + 157 037 | 0.9952 |
| nonanal | < 0.001 | < 0.001 | 0.01 | 0.01-2.5 | y = 1 x 10^7 x + 123 911 | 0.9942 |
| 2-nonenal | < 0.001 | < 0.001 | 0.05 | 0.05-2.5 | y = 1 x 10^7 x - 186 605 | 0.9925 |
| 3-methylbutanal | < 0.001 | 0.01 | 0.05 | 0.05-2.5 | y = 2 x 10^6 x + 1 x 10^6 | 0.9921 |
| benzaldehyde | 0.01 | 0.01 | 0.01 | 0.01-2.5 | y = 7 x 10^6 x + 17 925 | 0.9905 |
| 1-octen-3-ol | < 0.001 | 0.005 | 0.005 | 0.005-1.0 | y = 2 x 10^7 x + 227 470 | 0.9905 |
| 3-octen-2-one | < 0.001 | 0.005 | 0.005 | 0.005-2.0 | y = 2 x 10^7 x + 531 187 | 0.9950 |
| 2-pentylfuran | < 0.001 | 0.01 | 0.01 | 0.01-2.5 | y = 3 x 10^6 x - 87 598 | 0.9931 |
| 2,5-dimethylpyrazine | 0.005 | 0.05 | 0.05 | 0.05-2.0 | y = 1 x 10^6 x + 14 185 | 0.9942 |

In sample units (0.2 g in 2 mL): LOQ x 10 = ug/g isolate, i.e. hexanal LOQ 0.5 ug/g, 1-octen-3-ol
0.05 ug/g, 2-pentylfuran 0.1 ug/g; upper linear limit hexanal 25 ug/g, 1-octen-3-ol 10 ug/g.

### Table 3. "Profile of volatile compounds for the studied pea protein isolate at different pHs (mean of three repetitions +/- standard deviation)."

Unit: **mean peak area per g of sample** (A.U./g; full-scan TIC integration, not calibrated). Letters
= Tukey groups within a row; no SD = found in one of three repetitions; n.d. = not detected; LRI
on DB-WAX ("/" = below 3 min, not calculable). Bold in the original = semi-quantified compounds.
Literature column (which earlier papers reported the compound; X = none) omitted except where "X".
Sample = PPI1 (see §2).

| family | compound | CAS | RT min | LRI | pH 6.5 | pH 4.5 | pH 2 |
|---|---|---|---:|---:|---:|---:|---:|
| aldehyde | 3-methylbutanal | 590-86-3 | 2.4 | / | n.d. | 249 564 | 272 912 +/- 31 885 |
| | pentanal | 110-62-3 | 3.1 | / | 1 805 405 +/- 92 876 b | 3 969 920 +/- 410 505 a | 3 871 487 +/- 1 176 475 a |
| | **hexanal** | 66-25-1 | 5.3 | 1078 | **26 205 328 +/- 752 509 c** | **41 138 859 +/- 2 634 607 a** | **34 116 982 +/- 1 122 927 b** |
| | 2-hexenal | 505-57-7 | 9.6 | 1213 | 345 986 +/- 9 894 b | 736 319 +/- 88 963 a | 296 021 +/- 24 358 b |
| | heptanal | 111-71-7 | 8.6 | 1185 | 1 131 897 +/- 84 028 b | 2 546 586 +/- 355 942 a | 2 153 028 |
| | (E)-2-heptenal | 18829-55-5 | 14.0 | 1325 | 357 735 +/- 30 035 c | 2 006 121 +/- 247 892 a | 1 510 653 +/- 95 130 b |
| | octanal | 124-13-0 | 12.7 | 1293 | 1 121 496 +/- 70 638 b | 2 780 324 +/- 343 593 a | 990 280 +/- 860 436 b |
| | (E)-2-octenal | 2548-87-0 | 18.3 | 1432 | n.d. | 1 404 734 +/- 180 826 b | 2 175 206 +/- 72 577 a |
| | **nonanal** | 124-19-6 | 17.0 | 1400 | 5 465 053 +/- 412 073 ab | 8 175 080 +/- 1 137 161 a | 3 275 647 +/- 1 507 864 b |
| | **(E)-2-nonenal** | 18829-56-6 | 22.5 | 1539 | n.d. | 501 542 +/- 83 258 b | 680 233 +/- 66 018 a |
| | decanal | 112-31-2 | 21.3 | 1507 | n.d. | 425 076 +/- 77 308 | 274 347 |
| | (E)-2-decenal | 3913-81-3 | 26.3 | 1650 | n.d. | n.d. | 534 212 +/- 19 500 |
| | **benzaldehyde** | 100-52-7 | 21.8 | 1521 | 1 511 999 +/- 100 318 b | 2 264 672 +/- 331 285 a | 1 873 085 +/- 84 585 ab |
| alcohol | ethanol | 64-17-5 | 2.6 | / | 258 605 +/- 8 888 | n.d. | 332 070 +/- 84 559 |
| | 1-pentanol | 71-41-0 | 11.5 | 1262 | 502 134 +/- 60 348 b | 674 143 +/- 79 569 a | 562 929 +/- 38 614 ab |
| | **1-hexanol** | 111-27-3 | 15.8 | 1370 | **2 364 526 +/- 55 306 ab** | **2 654 065 +/- 204 139 a** | **2 112 372 +/- 90 927 b** |
| | 2-ethyl-1-hexanol | 104-76-7 | 21.4 | 1510 | n.d. | 1 425 349 +/- 132 803 b | 2 124 221 +/- 216 848 a |
| | 3-ethyl-3-hexanol (X) | 597-76-2 | 20.2 | 1479 | n.d. | 322 006 +/- 30 072 | n.d. |
| | 3,4-dimethyl-3-hexanol (X) | 19550-08-4 | 13.3 | 1308 | n.d. | 935 402 +/- 134 020 | n.d. |
| | 1-heptanol | 111-70-6 | 20.0 | 1475 | 248 472 +/- 9 704 a | n.d. | 385 807 +/- 105 604 a |
| | 1-octanol | 111-87-5 | 24.0 | 1579 | 466 763 +/- 35 326 a | 476 150 +/- 53 645 a | 584 164 +/- 51 390 a |
| | **1-octen-3-ol** | 3391-86-4 | 19.8 | 1470 | **1 031 787 +/- 65 657 b** | **1 453 934 +/- 149 255 a** | **1 134 591 +/- 130 594 b** |
| ketone | acetone | 67-64-1 | 1.7 | / | 341 103 +/- 9 898 a | 436 096 +/- 90 209 a | 460 636 +/- 183 006 a |
| | 2-heptanone | 110-43-0 | 8.5 | 1182 | 2 914 251 +/- 719 581 a | 3 384 973 +/- 1 279 843 a | 2 398 964 +/- 168 258 a |
| | 1-hepten-3-one (X) | 2918-13-0 | 13.0 | 1301 | n.d. | 861 709 +/- 178 335 | n.d. |
| | 3-methyl-2-heptanone | 2371-19-9 | 2.2 | / | n.d. | 327 462 | n.d. |
| | 6-methyl-5-hepten-2-one | 110-93-0 | 14.7 | 1343 | 282 010 +/- 13 561 | 293 682 | 275 776 +/- 21 363 |
| | 2-octanone | 111-13-7 | 12.5 | 1288 | 363 899 +/- 14 977 b | 393 851 +/- 39 817 b | 475 644 +/- 13 113 a |
| | 3-octanone | 106-68-3 | 11.2 | 1254 | n.d. | 262 761 | 296 389 |
| | **3-octen-2-one** | 18402-82-9 | 17.5 | 1412 | 609 220 +/- 174 395 a | 932 755 +/- 130 495 a | 692 621 +/- 103 697 a |
| | 3,5-octadien-2-one | 30086-02-3 | 23.9 | 1576 | 1 695 929 +/- 164 026 a | 1 353 879 +/- 228 926 ab | 1 083 916 +/- 31 733 b |
| | 2,3-octanedione | 585-25-1 | 14.4 | 1335 | 796 320 +/- 64 074 b | 1 100 100 +/- 111 882 a | 1 162 139 +/- 33 550 a |
| | 2-nonanone | 821-55-6 | 16.9 | 1397 | 715 439 +/- 71 852 a | 603 685 +/- 78 975 a | 611 501 +/- 13 248 a |
| | 2-decanone | 693-54-9 | 21.1 | 1502 | 476 536 +/- 76 890 a | 363 298 +/- 42 518 a | 381 605 +/- 3 068 a |
| | 1-decen-3-one (X) | 56606-79-2 | 13.0 | 1301 | n.d. | 937 292 | n.d. |
| furan | 2-ethylfuran | 3208-16-0 | 2.8 | / | 711 604 +/- 41 780 b | 860 919 +/- 99 991 b | UNREADABLE as printed ("2 030 27 +/- 88 2250" a; ~2.0 x 10^6) |
| | 2-n-butylfuran (X) | 4466-24-4 | 6.7 | 1127 | n.d. | 329 426 +/- 73 257 b | 628 836 +/- 5 984 a |
| | **2-pentylfuran** | 3777-69-3 | 10.4 | 1234 | **14 612 684 +/- 1 569 213 b** | **13 131 839 +/- 1 809 401 b** | **21 303 106 +/- 2 213 382 a** |
| | cis/trans-2-(2-pentenyl)furan (X) | 70424-13-4 | 13.2 | 1306 | n.d. | n.d. | 739 349 +/- 235 952 |

39 compounds: 13 aldehydes, 9 alcohols, 13 ketones, 4 furans (counts reproduce). The six
compounds not previously reported (X) were "retrieved only at acidic pH extraction". No pyrazine of
any kind was detected.

### Table 4. "Semi-quantification of the compounds of interest in two pea protein isolates at different pH (mean of three repetition +/- standard deviation, in ug of compound/g of sample)."

n.d. = not detected; n.q. = not quantified (below LOQ); > uloq = above the upper limit of
quantification. Letters = Tukey groups within a row and within one isolate.

| compound | PPI1 pH 6.5 | PPI1 pH 4.5 | PPI1 pH 2.0 | PPI2 pH 6.5 | PPI2 pH 4.5 | PPI2 pH 2.0 |
|---|---:|---:|---:|---:|---:|---:|
| **hexanal** | **5.1 +/- 0.1 c** | **8.1 +/- 0.5 a** | **6.7 +/- 0.2 b** | **3.4 +/- 0.1 c** | **5.6 +/- 0.4 b** | **7.3 +/- 0.4 a** |
| nonanal | 0.97 +/- 0.08 ab | 1.5 +/- 0.2 a | 0.6 +/- 0.3 b | 0.67 +/- 0.07 a | 0.8 +/- 0.2 a | 0.6 +/- 0.4 a |
| (E)-2-nonenal | n.d. | 0.26 +/- 0.06 a | 0.313 +/- 0.008 a | n.d. | 0.28 +/- 0.02 b | 0.56 +/- 0.09 a |
| 3-methylbutanal | n.d. | n.q. | n.q. | 4.6 +/- 0.4 a | 5.7 +/- 0.8 a | 2 +/- 1 b |
| benzaldehyde | 0.41 +/- 0.03 b | 0.62 +/- 0.09 a | 0.51 +/- 0.02 ab | 1.6 +/- 0.1 b | 1.5 +/- 0.1 b | 2.8 +/- 0.3 a |
| **1-octen-3-ol** | **0.07 +/- 0.01 b** | **0.12 +/- 0.01 a** | **0.08 +/- 0.01 b** | **0.20 +/- 0.02 a** | **0.29 +/- 0.06 a** | **0.27 +/- 0.03 a** |
| 3-octen-2-one | n.q. | n.q. | n.q. | 0.21 +/- 0.02 b | 0.39 +/- 0.09 b | 0.7 +/- 0.2 a |
| **2-pentylfuran** | **10 +/- 1 b** | **9 +/- 1 b** | **14 +/- 1 a** | **12.0 +/- 0.8 a** | **11.5 +/- 0.7 a** | **> uloq** |
| 2,5-dimethylpyrazine | n.d. | n.d. | n.d. | n.d. | n.d. | n.d. |

Other numbers in the text: "hexanal release was found 59 % higher with extraction using pH 4.5
than with pH 6.5" (8.1/5.1 = 1.59, PPI1); "5.1 ug of hexanal/g of sample were initially free in the
matrix and 3.0 ug/g of sample were bound to the matrix" (the authors' free/bound split, PPI1);
Heng et al. 2004 cited: vicilin binds 1-17 % of ketones vs 75-88 % of aldehydes; legumin binds no
ketones.

## 4. Numbers the repository can use

Molar conversions: hexanal MW 100.16, 1-hexanol 102.17, 2-pentylfuran 138.21, 1-octen-3-ol 128.21.
"Per g isolate" is per g of powder as received (85 % protein on dry matter; moisture not given).

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| **Hexanal, free (neutral extraction), PPI1 / PPI2** | 5.1 +/- 0.1 / 3.4 +/- 0.1 (= 51 / 34 nmol/g) | ug per g isolate | 10 % w/v in water, pH 6.5, 40 C, 30 + 60 min HS-SPME, external calibration in water | Table 4 | level_only (the carried-over hexanal of a commercial pea isolate) |
| Hexanal, maximum released (any pH), PPI1 / PPI2 | 8.1 +/- 0.5 (pH 4.5) / 7.3 +/- 0.4 (pH 2.0) | ug/g | same | Table 4 | level_only |
| Hexanal pH 4.5 / pH 6.5 | 1.59 (PPI1), 1.65 (PPI2) | x | same lot, same fibre | Table 4 | within_study_ratio (the authors' "59 %") |
| Hexanal pH 2.0 / pH 6.5 | 1.31 (PPI1), 2.15 (PPI2) | x | | Table 4 | within_study_ratio (direction differs between lots) |
| Hexanal bound fraction at pH 6.5 (authors' split) | 3.0 of 8.1 = 37 % bound, 63 % free (PPI1) | % | | §3.3.3 | within_study_ratio (assumes pH 4.5 releases everything) |
| **1-Octen-3-ol, PPI1 / PPI2** | 0.07 +/- 0.01 / 0.20 +/- 0.02 (pH 6.5); 0.12 / 0.29 (pH 4.5); 0.08 / 0.27 (pH 2.0) | ug/g | same | Table 4 | level_only; PPI1 values sit at 1-2x the LOQ (0.05 ug/g) and the printed Table 2 equation does not reproduce them from Table 3 (flag 5) |
| **2-Pentylfuran, PPI1 / PPI2** | 10 +/- 1 / 12.0 +/- 0.8 (pH 6.5); 9 / 11.5 (pH 4.5); 14 / > 25 (pH 2.0) | ug/g | same; > uloq = above 2.5 ppm assay = > 25 ug/g | Table 4 | level_only; **2-pentylfuran is the largest beany volatile by mass in this isolate, 2-3x hexanal** |
| 2-Pentylfuran pH 2.0 / 6.5 | 1.4 (PPI1); > 2.1 (PPI2) | x | | Table 4 | within_study_ratio (no release at pH 4.5, release at pH 2) |
| **1-Hexanol** | NOT semi-quantified; peak area/g 2 364 526 / 2 654 065 / 2 112 372 at pH 6.5 / 4.5 / 2 (PPI1) | A.U./g | | Table 3 | level_only in area units — not convertible to ug/g (no hexanol calibration) |
| 1-Hexanol : hexanal peak-area ratio, pH 6.5 | 0.090 | area/area | full-scan TIC, uncorrected response factors | Table 3 | within_study_ratio in area units only; contrast raw pea milk (Zhang 2020b) where hexanol exceeds hexanal 2.3:1 by mass — in the finished isolate hexanal dominates (flag 6) |
| 1-Hexanol pH 4.5 / 6.5 (area) | 1.12 | x | | Table 3 | within_study_ratio |
| Nonanal, PPI1 / PPI2 | 0.97 / 0.67 (pH 6.5); 1.5 / 0.8 (pH 4.5); 0.6 / 0.6 (pH 2.0) | ug/g | | Table 4 | level_only |
| (E)-2-Nonenal | n.d. at pH 6.5; 0.26 / 0.28 (pH 4.5); 0.313 / 0.56 (pH 2.0) | ug/g | | Table 4 | level_only; a fully bound aldehyde at neutral pH |
| Benzaldehyde, PPI1 / PPI2 | 0.41 / 1.6 (pH 6.5); 0.62 / 1.5; 0.51 / 2.8 | ug/g | | Table 4 | level_only |
| 3-Methylbutanal | n.d. (PPI1); 4.6 / 5.7 / 2 (PPI2) | ug/g | | Table 4 | level_only; Strecker aldehyde present in one lot only — the authors read it as lot age |
| 3-Octen-2-one | n.q. (PPI1); 0.21 / 0.39 / 0.7 (PPI2) | ug/g | | Table 4 | level_only |
| 2,5-Dimethylpyrazine | n.d. in both lots at all pH (< 0.05 ppm assay = < 0.5 ug/g) | ug/g | | Table 4 | level_only (null) — an unheated pea isolate carries no pyrazine above 0.5 ug/g |
| Sum of nine semi-quantified beany volatiles, pH 6.5 | ~16.6 (PPI1); ~18.1 (PPI2), 2-pentylfuran + hexanal = 91-85 % of it | ug/g | | derived from Table 4 | level_only |
| Charge arithmetic for a recipe (not a fit) | at 30 g protein/L (35 g isolate/L): hexanal 120-180 ug/L free, 2-pentylfuran 350-420 ug/L, 1-octen-3-ol 2.5-7 ug/L | ug/L | 35 g/L x Table 4 pH 6.5 values | derived | derived level — same order as raw pea milk hexanal 164 ug/L (Zhang 2020b) |
| Isolate identity | Roquette spray-dried PPI, 85 % protein (dry matter), mainly globulins | — | | §2.4.1 | protocol |
| Method precision | ~15 % single-fibre CV; 6-68 % between fibres; LOQ hexanal 0.05 ppm assay (0.5 ug/g) | — | | abstract, Tables 1-2 | method |
| LOX activity, lipid content, time course | NOT MEASURED | — | | — | — |

## 5. Flags

1. **Semi-quantification by external calibration in pure water.** No internal standard, no
   matrix matching, no recovery correction: the isolate suspension binds aldehydes (the paper's
   own point), so the headspace partition in the sample is lower than in the calibration water
   and every ug/g is an under-estimate of total content by an unknown, compound- and pH-dependent
   factor. Levels are lower bounds on the powder's content; within-lot pH ratios are sound.
2. **Two lots, no provenance.** PPI1 and PPI2 differ 1.5x in hexanal and qualitatively in
   3-methylbutanal (n.d. vs 4.6 ug/g) and 3-octen-2-one; the authors guess age. Any repository
   row should carry both lots as a band, not a mean.
3. **Residual LOX during the 90-min, 40 C headspace step is uncontrolled.** Gao 2020 shows a
   finished pea isolate retains LOX; this isolate's lipid content is not given. At pH 6.5 in water
   the vial is a mild incubation; at pH 2 LOX is inactive. If some hexanal formed in the vial, the
   "free at pH 6.5" number is inflated and the pH 4.5 / 6.5 ratio deflated. The paper offers no
   heat-killed control.
4. **Extraction pH is a release pH, not a process pH.** The isolate itself was made at Roquette's
   (unstated) conditions; the pH series probes protein-volatile binding in the vial. It says
   nothing about LOX kinetics and must not be used to fit any rate; it can inform the
   `protein_matrices.yml` binding side (aldehydes ~37 % bound at pH 6.5 by the authors' split).
5. **1-Octen-3-ol arithmetic does not close.** Table 3 area/g x 0.2 g through the printed Table 2
   equation gives negative or ~0.03 ug/g where Table 4 prints 0.07-0.12; hexanal, nonanal,
   benzaldehyde and 2-pentylfuran close to the printed precision. The printed 1-octen-3-ol slope /
   intercept are rounded too coarsely near its LOQ. Treat PPI1 1-octen-3-ol as order-of-magnitude
   (0.05-0.15 ug/g).
6. **1-Hexanol is area-only.** No calibration was run for it, so the hexanol : hexanal contrast
   with raw pea milk (Zhang 2020b, 2.3:1 by mass) is a peak-area statement (0.09:1) with
   uncorrected, different response factors; the direction (hexanal >> hexanol in the isolate,
   hexanol > hexanal in raw milk) is probably real but the magnitude is not established.
7. **One Table 3 cell is unreadable** (2-ethylfuran, pH 2; garbled in the PDF text layer itself);
   2-ethylfuran is not a repository species, no action.
8. **Table 3 = PPI1 is an inference** verified by arithmetic (§2), not a caption statement.
9. **Table 1 "CV between two fibres" for hexanal (18.85 %)** does not follow from the two printed
   means (4 373 173 vs 4 207 065, ~4 % apart); the column's definition is not given. Not
   load-bearing.
10. **Units and dilution:** 0.2 g in 2 mL; assay ppm x 10 = ug/g isolate; upper linear limits
    (hexanal 2.5 ppm = 25 ug/g) were not approached except 2-pentylfuran in PPI2 at pH 2 (> uloq).
