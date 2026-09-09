# Cerny & Davidek 2004 — EXTRACTION (alpha-mercaptoketone formation from cysteine + [1-13C]ribose, pH 5, 95 C / 4 h)
### A labelling-position paper: it fixes where ribose C-1 ends up in MFT, FFT and the mercaptoketones. It contains NO yields or concentrations.

**Source on disk:** `data/articles/cerny2004.pdf` (owner's download, 2026-09-07). Read-only extraction from
`pdftotext -layout`; the single table (MS spectra + labelling distribution) has a clean text layer. Repo
status before this dossier: FIT_HOLDOUT_DECLARATION lists "Cerny 2004 (in-situ branching 54:46 etc.) =
FIT"; k3 §A.10 carries the 54:46 butanedione split, 65:35 thiazole and 87:13 methylthio splits, and the
"C2+C3 route to MFT was not relevant at 95 C" claim. This dossier confirms those numbers against the PDF
and records that the paper holds nothing else quantitative.

## 0. Identity

| field | value |
|---|---|
| Title | "alpha-Mercaptoketone Formation during the Maillard Reaction of Cysteine and [1-13C]Ribose" |
| Authors | Christoph Cerny and Tomas Davidek (Nestle Research Center, Lausanne) |
| Venue | J. Agric. Food Chem. 2004, 52, 958-961 |
| DOI | 10.1021/jf035265m |
| Companion | ref 20 = Cerny & Davidek 2003, JAFC 51, 2714-2721 ([13C5]ribose; the repo's cerny2003 HOLD-OUT) |

## 1. Why it matters

The task brief hoped for mercaptoketone yields. **There are none.** The paper is purely a positional
13C-labelling study: which ribose carbon becomes which product carbon, and in what proportion, for 12
volatiles at cooking conditions (95 C, 4 h, pH 5). What it gives the model is topology: (i) 3-mercapto-
2-pentanone carries the intact ribose chain with ribose C-1 as its C-1 (the 1,4-dideoxypento-2,3-diulose
route), (ii) 2-mercapto-3-pentanone is NOT formed at 95 C, (iii) butane-2,3-dione (hence 3-mercapto-
2-butanone) arises 54 : 46 from loss of ribose C-1 vs C-5, i.e. two parallel C4 routes, (iv) MFT carries
ribose C-1 as its 2-methyl carbon, (v) FFT carries ribose C-1 as its CH2SH carbon (via 3-deoxypentosone
and furfural), and (vi) the hydroxyacetaldehyde + mercapto-2-propanone (C2 + C3) route to MFT "was not
relevant under the reaction conditions used in this study".

## 2. Methods as they matter to a model

- **Charges:** cysteine 5.84 mg (MW 121.16 = 48.2 µmol) + ribose 21.75 mg (MW 150.13 = 144.9 µmol) in
  **472 mg of 0.5 mol/L potassium phosphate buffer, pH 5.00** (~0.47 mL) -> roughly **100 mmol/L cysteine,
  300 mmol/L ribose, 1 : 3** (the authors state only the 1 : 3 ratio). Second run identical with
  [1-13C]ribose (99 % enrichment).
- **Vessel / heating:** 2 mL glass vials, heated metal block (Reacti-Therm), **95 C, 4 h**. Stirring not
  stated for this paper (the block is a stirring/heating module). Atmosphere: sealed vial headspace.
- **Replicates:** one labelled and one unlabelled run; no replication stated.
- **Analysis:** HS-SPME directly from the reaction vial, PDMS-DVB 65 µm, 60 min at 20 C, no agitation;
  GC 6890A, HP-PONA 50 m x 0.20 mm x 0.50 µm, 35 -> 240 C at 6 C/min; MSD 5973, EI 70 eV, m/z 28-350.
  RI on OV-1.
- **Quantification: NONE.** No internal standard, no peak areas, no concentrations. The only numbers are
  isotopomer ratios from integrated molecular-ion (and fragment-ion) signals: "Integration of the
  molecular ion signals at m/z 104 and 105 indicates 54% unlabeled (5a) and 46% monolabeled
  3-mercaptobutan-2-one"; "The ratio of the corresponding integrated peaks in the ion chromatograms was
  found to be 74:26 (m/z 43:44)"; thiazole 65 % labelled; 2-methyl-3-(methylthio)furan 87 : 13.

## 3. The table re-typed

### Table 1. "MS Spectra of Compounds Formed from Cysteine and [1-13C]Ribose"

Columns in the paper: no.; compound; RI (OV-1); reaction cysteine + ribose (reference) m/z (%); reaction
cysteine + [1-13C]ribose m/z (%); 13C-labeling distribution (no / 13C1 / 13C2). The two m/z lists are
spectral data (re-typed only for the model-relevant species); the last column is the quantitative
content.

| no. | compound | RI (OV-1) | labelling distribution (unlabelled / 13C1 / 13C2), % | position of the 13C (text) |
|---:|---|---:|---|---|
| 1 | furan | 500 | 100 / 0 / 0 | ribose C-1 lost |
| 2 | 2-methylfuran | 586 | 0 / 100 / 0 | methyl carbon |
| 3 | thiazole | 709 | 35 / 65 / 0 | C-2 (from [13C]formaldehyde = ribose C-1); 35 % from unlabelled formaldehyde |
| 4 | 2-methylthiophene | 756 | 0 / 100 / 0 | methyl carbon |
| 5 | **3-mercaptobutan-2-one** | 782 | **54 / 46 / 0** | of the 46 % labelled: 26 % at C-1 (5b), 20 % at C-4 (5c), from the m/z 43 : 44 = 74 : 26 fragment ratio |
| 6 | **furan-2-carbaldehyde (furfural)** | 801 | 0 / 100 / 0 | aldehyde carbon (from 3-deoxypento-1,2-diulose) |
| 7 | **2-methyl-3-furanthiol (MFT)** | 850 | 0 / 100 / 0 | 2-methyl carbon; C-5 of the ring unlabelled (loss of H12CO, m/z 29) |
| 8 | **3-mercaptopentan-2-one** | 871 | 0 / 100 / 0 | C-1 (base peak m/z 44 = 13CH3CO+; M+ 119) — intact ribose chain via 1,4-dideoxypento-2,3-diulose |
| 9 | **2-furfurylthiol (FFT)** | 883 | 0 / 100 / 0 | CH2SH carbon (M+ 115 -> 82 -> 54: loss of SH then 12CO from C-5) |
| 10 | 2-methyl-3-(methylthio)furan | 933 | 0 / 87 / 13 | methyl always labelled; methylthio carbon labelled in 13 % (so only part of the ribose-derived S-CH3 comes from C-1) |
| 11 | 3-thiophenethiol | 940 | 100 / 0 / 0 | from cysteine, not ribose |
| 12 | 2-methyl-3-thiophenethiol | 1042 | 0 / 100 / 0 | methyl carbon |

Reference-spectrum m/z lists for the species the model carries (unlabelled run / labelled run):
- 5: 43 (100), 58 (46), 60 (29), 61 (82), 104 (48) / 43 (100), 44 (26), 58 (69), 60 (36), 61 (70), 62 (30), 104 (33), 105 (34)
- 6: 29 (8), 39 (34), 67 (6), 95 (94), 96 (100) / 30 (5), 38 (11), 39 (33), 67 (11), 96 (98), 97 (100)
- 7: 43 (13), 71 (18), 85 (25), 113 (24), 114 (100) / 44 (15), 71 (16), 86 (25), 114 (25), 115 (100)
- 8: 39 (13), 41 (60), 43 (100), 47 (40), 74 (59), 75 (71), 118 (25) / 39 (13), 41 (59), 44 (100), 47 (41), 74 (58), 75 (68), 119 (25)
- 9: 53 (39), 81 (100), 114 (35) / 54 (36), 82 (100), 115 (33)

**Species NOT detected:** 2-mercapto-3-pentanone ("the isomer 2-mercapto-3-pentanone was not detected",
consistent with ref 20 at the same conditions). No 4-hydroxy-5-methyl-3(2H)-furanone, no bis(2-methyl-
3-furyl) disulfide, no H2S in the list (SPME/GC window).

## 4. What the repo could take

No FIT rows in the yield sense — there is not a single concentration or peak area in the paper. What it
gives are step-level branching fractions and route identities at **95 C / 4 h / pH 5.00 / 0.5 M phosphate
/ cys : rib = 1 : 3 (~100 : 300 mM)**:

| constraint | value | what it pins |
|---|---|---|
| butane-2,3-dione (-> 3-mercapto-2-butanone) from loss of ribose C-1 vs C-5 | **54 : 46** (M+ 104 : 105) | two parallel C4 routes: retro-Mannich of the cysteine Amadori/5-deoxyosone (loses C-1, Fig. 2) vs retro-aldol of 1,4-dideoxypento-2,3-diulose or 1-deoxypentosone (loses C-5, Fig. 3). Already FIT in the repo. |
| within the labelled 3-mercaptobutan-2-one, label at C-1 vs C-4 | **26 : 20** (74 : 26 on m/z 43 : 44) | the two labelled isomers 5b / 5c; the diketone + H2S step does not discriminate |
| 3-mercapto-2-pentanone | 100 % singly labelled, label at C-1 | intact C5 chain, 1,4-dideoxypento-2,3-diulose route; NF is not needed |
| 2-mercapto-3-pentanone | **not detected** at 95 C | the NF-diagnostic isomer is a high-temperature product (Whitfield 1999: 77.5 µg/10 mg NF at 140 C; Cerny 2007: 0 at pH <= 5.5 even at 145 C, 19 at pH 6) |
| MFT | 100 % singly labelled at the 2-methyl carbon | ribose C-1 -> methyl; **C2 + C3 recombination "was not relevant"** at 95 C. Combined with Cerny 2003 (intact C5, 49 : 46 with 13C5), the 95 C MFT is an intact-skeleton, C-1-methyl product |
| FFT | 100 % singly labelled at CH2SH | ribose C-1 -> furfural CHO -> FFT CH2SH: the furfural + H2S step is confirmed as the FFT route at 95 C |
| furfural | 100 % labelled at CHO | 3-deoxypentosone route, no fragmentation |
| thiazole | 65 % labelled at C-2 | 65 : 35 formaldehyde from ribose C-1 vs elsewhere |
| 2-methyl-3-(methylthio)furan | 87 : 13 (methylthio C unlabelled : labelled) | only part of the ribose-derived S-methyl comes from C-1 |
| furan | 0 % labelled | C-1 lost |
| 3-thiophenethiol | 0 % labelled | cysteine-derived, not ribose |

Directional claim for the temperature axis (text): "Recently, Schieberle and co-workers showed that at a
higher reaction temperature (180 °C) and low water content (10%) butane-2,3-dione is not formed from an
intact carbohydrate chain, but mainly from C3/C1 recombination (21). Under cooking conditions (95 °C, 4 h),
as used in the present study, fragmentation was negligible." So the C4 diketone's origin (and, by the
authors' reading, MFT's) switches from intact-chain at 95 C to recombination at 180 C dry — the repo's
"C2+C3 lane is temperature-scoped" owner call rests on this sentence plus Hofmann 1998 T10.

## 5. Caveats

1. **Zero quantitative yield information.** Any "mercaptoketone yield" attributed to this paper elsewhere
   in the repo would be a mis-citation; the yields the repo uses for mercaptoketones come from Whitfield
   1999 (140 C) and Hofmann 1998 (145 C).
2. **Single run each**, no replication, isotopomer ratios from integrated ion signals without stated
   uncertainty; treat 54 : 46 as 50 +- 5.
3. The 26 : 20 split inside the labelled 3-mercaptobutan-2-one rests on a fragment-ion ratio (m/z 43 : 44)
   and on the assumption (from ref 20) that the C4 chain does not fragment; it is a derived number.
4. **95 C / 4 h**, 50 C below the fit panel; the paper itself says the routes change with temperature.
   The branching fractions are FIT only at Cerny's conditions, as the existing declaration already notes
   for cerny2003.
5. The buffer is given by mass (472 mg), not volume; the ~100 / 300 mM concentrations above are the
   dossier's arithmetic (density ~1), twice the Cerny 2003 loading if that paper used 1 mL — check before
   using concentration in any comparison.
6. HS-SPME at 20 C on PDMS-DVB, no agitation: fine for isotopomer ratios (same compound, same partition),
   not for anything between compounds.
