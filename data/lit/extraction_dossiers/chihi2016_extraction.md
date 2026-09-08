# Chihi, Mession, Sok & Saurel 2016 — EXTRACTION (pea globulins and BLG, free SH and S-S by Ellman's in 1.5 M GdnCl, native and 85 C / 60 min)
### The only paper on disk that titrates pea globulin AND beta-lactoglobulin side by side under one protocol; Table 1 is the whole payload.

**Source on disk:** `data/articles/chihi2016.pdf` (owner's download, 2026-09-08). Read from the pdftotext
layer in the scratchpad (`articles/chihi2016.txt`); the text layer of Table 1 is clean and is re-typed in
full below, as are Tables 2 and 3 in condensed form. This is the ASAP / "XXXX" version (no volume or
page numbers printed). Repo status before this dossier: `data/species/protein_matrices.yml` carries BLG
with computed 0.05446 mmol/g free thiol and 0.10892 mmol/g disulfide; no pea entry.

## 0. Identity

| field | value |
|---|---|
| Title | "Heat-Induced Soluble Protein Aggregates from Mixed Pea Globulins and β-Lactoglobulin" |
| Authors | Mohamed-Lazhar Chihi, Jean-luc Mession, Nicolas Sok, Rémi Saurel (UMR PAM, AgroSup Dijon; Chihi also Alger) |
| Venue | J. Agric. Food Chem., ASAP 2016 (received 7 Jan, accepted 20 Mar 2016); volume/pages not printed |
| DOI | 10.1021/acs.jafc.6b00087 |
| Naming | "Glob" = laboratory pea globulin isolate (7S vicilin/convicilin + 11S legumin); "βlg" = laboratory-purified beta-lactoglobulin from WPI (variant not stated); "S− free" = free sulfhydryl; "S−S" = disulfide; H0 = ANS surface-hydrophobicity slope x 10^6. Mixtures are βlg/Glob weight ratios. |
| Companion | Mession et al. 2013, JAFC 61, 1196-1204 (pea globulin heat aggregation; source of the heating protocol and the "comparable" native Glob SH values); not on disk |

## 1. Why it matters

The matrix layer has BLG's sites from its sequence (1 free Cys, 2 S-S per 18 362 Da) and nothing for
pea. This paper measures, with one Ellman protocol, (a) BLG free SH 42.5 and S-S 102.5 µmol/g protein,
a direct check on the computed 54.5 and 108.9, and (b) pea globulin free SH 2.1 and S-S 4.2 µmol/g
protein, the first pea site density on disk, together with both after 85 C / 60 min at 2 wt % and low
ionic strength. It also states the composition of the pea fraction (legumin 35 %, vicilin 45 %,
convicilin 9 %, lipoxygenase 7 %), which matters because vicilin carries no cysteine and legumin
carries nearly all of it.

## 2. Methods as they matter to a model

- **Pea globulin isolate (lab):** smooth yellow pea flour (Roquette), defatted (petroleum ether +
  ethanol, 4 C, 1 h); albumins removed by 0.1 M acetate pH 4.9 (2 h, 20 C, 12 000 g); globulins
  extracted from the washed pellet with 0.1 M Na2HPO4 + 5 % (w/v) K2SO4, pH 8, 1:10 (w/v), 2 h, 4 C;
  UF/DF (10 kDa) against 5 mM ammonium carbonate pH 7.2; freeze-dried, -20 C. DSC: Td 76 C, ΔHd 14.4
  J/g protein (native). Composition: "protein and fat contents were around 94 and <1 wt %" (dry basis);
  Kjeldahl N x 6.25. SDS-PAGE: legumin ~35 %, vicilin ~45 %, convicilin ~9 %, lipoxygenase ~7 % of band
  intensity.
- **BLG isolate (lab):** Promilk 802 FB WPI (Ingredia) 10 wt %, TCA to pH 2 (30 min, 4 C), 12 000 g,
  dialysis (10 kDa) 3x water then 3x 5 mM NaCl, 20 000 g, pH 7.2, freeze-dried. DSC Td ~71 C, ΔHd 7.9
  J/g protein. "∼92 wt % protein, ∼3 wt % ashes, and ∼0.8 wt % fats on a dry basis"; Kjeldahl N x 6.38.
  SDS-PAGE: monomer band + "negligible traces of non-native dimeric βlg". Genetic variant NOT stated
  (commercial WPI, so an A/B mixture is the default expectation).
- **Solutions and heating:** stock "2 wt %" in 10 mM sodium phosphate + 5 mM NaCl, pH 7.2, 0.02 %
  azide, 24 h 4 C, 12 000 g, 0.45 µm; mixtures 0/100, 30/70, 50/50, 70/30, 100/0 βlg/Glob at 2 wt %
  total, 2 h 25 C. Sealed tubes, 40 -> 85 C at 1 C/min, 85 C for 60 min, ice 10 min. Insolubles removed
  12 000 g 20 min; soluble protein after heating ~98 % of total for all ratios. DSC after heating: no
  endotherm.
- **SH method, verbatim:** "The free sulfhydryl (S−) and disulfide bond (S−S) contents of the protein
  samples were measured using 5,5′-dithiobis(2-nitrobenzoic acid) (DTNB), namely, Ellman's reagent. For
  S−free content determination, protein solutions at 2 wt % concentration were extensively dialyzed
  against a 0.1 M phosphate buffer at 4 °C, pH 7.5 (ratio protein solution-to-buffer 1:20). The reaction
  mixture in the phosphate buffer was prepared by mixing 650 μL of the protein sample with 750 μL of
  guanidium hydrochloride (3 M GdnCl) and 100 μL of DTNB solutions (1 mM), stirred vigorously, and kept
  in the dark for 10 min at 25 °C. Absorbance was measured at 412 nm using a molar extinction coefficient
  ε for DTNB of 12900 M−1/cm, as calculated from calibration curves with N-acetylcysteine (NAC) in the
  range of 0−60 μM. For S−S estimation, total sulfhydryl groups S−total (S−free + S−S) of the different
  protein samples were recovered by treatment with 20 mM DTT added as powder, for 2 h at 25 °C. Further
  extensive dialysis with degassed phosphate buffer allowed the elimination of the excess of DTT, with
  limited risk of reoxidation of the released S−. The S−S content was calculated as follows (eq 1):
  [S−S] = ([S−]total − [S−]free)/2. The free S− and S−S contents were expressed as μmol S−/g protein and
  μmol S−S/g protein, respectively."
  So: DTNB, final GdnCl 1.5 M (750 µL of 3 M in 1.5 mL), pH 7.5, 10 min, 25 C; "free" = denaturant-
  exposed free SH (not native-surface SH); total after 20 mM DTT then dialysis; S-S by the standard
  half-difference. ε = 12 900 M^-1 cm^-1 (calibrated with NAC; Shimada used 13 600).
- **Reagent stoichiometry as printed (flag, see §5.1):** 100 µL x 1 mM DTNB = 0.10 µmol per 1.5 mL
  (67 µM). 650 µL of a 2 wt % BLG solution holds ~13 mg protein x 42.5 µmol/g = 0.55 µmol SH (free) or
  x 247 µmol/g = 3.2 µmol SH (after DTT). The assay as written is therefore DTNB-limited for BLG by 5x
  (free) and 30x (total) unless the samples were diluted into the 0-60 µM NAC calibration range; a
  dilution step is implied by the calibration range but not stated.
- **Surface hydrophobicity:** ANS, 0.004-0.02 wt % protein in 10 mM Na2HPO4 pH 7.2, ex 390 / em 470 nm,
  H0 = slope of intensity vs concentration, n = 3.
- **Replicates:** Table 1 "mean ± standard deviation, calculated from three repetitions".
- **Units in the repo:** all Table 1 SH values are per g PROTEIN (Kjeldahl, 6.25 Glob / 6.38 BLG);
  1 µmol/g = 0.001 mmol/g. Powder basis: x 0.94 (Glob), x 0.92 (BLG), dry basis.

## 3. Tables re-typed

### Table 1. "Surface Hydrophobicity (H0) and Sulfhydryl Group (S−) and Disulfide Bridge (S−S) Contents of (a) Unheated and (b) Heated (85 °C for 60 min) Single-Protein Samples and βlg/Glob Mixtures at Different Weight Ratios"

Footnote a: "All protein samples were prepared at 2 wt % total protein samples, in 5 mM NaCl at pH 7.2.
All results are given as the mean ± standard deviation, calculated from three repetitions. Mean values
bearing the same letter (a−e) are not significantly different (p > 0.05). nd, not determined."
Footnotes b, c: methods as above.

| sample | H0 (slope x 10^6) | free S− (µmol/g protein) | S−S (µmol/g protein) |
|---|---:|---:|---:|
| **(a) Unheated** | | | |
| Glob | 2.4 ± 0.1 | **2.1 ± 0.1** | **4.2 ± 0.1** |
| 30/70 | 2 ± 0.1 | nd | nd |
| 50/50 | 1.7 ± 0 | nd | nd |
| 70/30 | 1.3 ± 0.1 | nd | nd |
| βlg | 1.1 ± 0.1 | **42.5 ± 0.2** | **102.5 ± 0.3** |
| **(b) Heated, 85 C / 60 min** | | | |
| Glob | 3.1 ± 0.1 a | **7.1 ± 0.2 a** | **1.8 ± 0.3 a** |
| 30/70 | 3.3 ± 0 b | 10.7 ± 0.3 ab | 35.8 ± 0.1 b |
| 50/50 | 3.4 ± 0.1 bc | 12 ± 0.4 b | 59 ± 0.3 c |
| 70/30 | 3.4 ± 0.1 bc | 14 ± 0.2 bc | 89 ± 0.1 d |
| βlg | 3.6 ± 0.1 c | **17.2 ± 0.3 c** | **114.1 ± 0.1 e** |

Running-text restatements: βlg free S− "from 42 to 17", S−S "from 102 to 114"; Glob: ">2-fold decrease
in S−S content, whereas the S− content increased markedly"; Liu et al. (ref 45) cited for ~28 µmol/g
accessible thiol in BLG aggregates after 85 C / 15 min at pH 7.

Half-cystine bookkeeping (this dossier's arithmetic, free + 2 x S-S):
- βlg: 42.5 + 205.0 = **247.5** native; 17.2 + 228.2 = **245.4** heated (conserved to 1 %).
- Glob: 2.1 + 8.4 = **10.5** native; 7.1 + 3.6 = **10.7** heated (conserved to 2 %).
- Mixtures, heated, weighted from the single-protein heated values: 30/70 expected free 0.3 x 17.2 + 0.7
  x 7.1 = 10.1 (printed 10.7); S-S 0.3 x 114.1 + 0.7 x 1.8 = 35.5 (35.8). 50/50: 12.2 (12); 58.0 (59).
  70/30: 14.2 (14); 80.4 (**89**). The 70/30 S-S exceeds additivity by ~9 µmol/g; the others are additive
  within 1 µmol/g.

### Table 2. Hydrodynamic diameter Dh (nm), DLS at 25 C, 0.1 wt % in water (condensed)

| sample | unheated | heated 85 C / 60 min |
|---|---|---|
| Glob | 16.9 ± 5.6 (56 ± 1.5 %) and 94.5 ± 1.5 (43.9 ± 0.4 %) | 69.3 ± 5.6 (63.1 %) and 151.9 ± 9.8 (36.9 %) |
| 30/70 | nd | 43.2 ± 0.2 (28.9 %) and 111 ± 10.1 (71.1 %) |
| 50/50 | nd | 37 ± 3.6 (15.6 %) and 95 ± 8.32 (84.4 %) |
| 70/30 | nd | 27 ± 0.3 (2.3 %) and 94.9 ± 0.1 (97.7 %) |
| βlg | 6 ± 0.4 (100 %) | 38 ± 0.07 (100 %) |

### Table 3. SEC-HPLC peak areas (280 nm), condensed
Native Glob: G1 660-443 kDa 23 %, G2 443-43 kDa 35 %, G3 43-29 kDa 17 %, G4 29-13.7 kDa 14 %, G5 < 13.7
kDa 10 %. Native βlg: one peak (45-14 kDa) 100 %. Heated Glob: G'1 (> 600 kDa) 11.9 %, G'2 (600-75)
49.5 %, G'3 (75-29) 30.4 %, G'4 (< 13.7) 8.2 %. Heated βlg: β'1 (> 2000 kDa) 100 %. Heated mixtures
30/70 / 50/50 / 70/30: > 2000 kDa 47.5 / 66.5 / 78.2 %; 2000-200 kDa 32.3 / 20.1 / 13.3 %; 75-29 kDa
20.2 / 13.4 / 8.5 %. ~50 % of heated single-Glob aggregates were lost on the 0.45 µm pre-filter (text).

### SDS-PAGE densitometry (text)
Heated single βlg: ~70 % of protein in > 200 kDa disulfide-bonded aggregates. Heated single Glob: ~10 %
of legumin subunits (Lα +11 %, Lβ1 +10 %, Lβ2 +11 % between NR and R lanes) disulfide-bonded; vicilin
patterns identical NR vs R (no covalent participation). Mixtures: non-migrating disulfide-bonded
aggregates ~33 / 37 / 55 % of polypeptides (30/70, 50/50, 70/30); legumin participation 70-95 %;
convicilin 20-30 % (covalent, only in the presence of βlg).

## 4. Site densities the repository can use

Protein basis: per g protein (Kjeldahl) as printed; no conversion assumption.

| matrix | quantity | value ± sd | unit as printed | mmol per g PROTEIN | conditions | source | evidence |
|---|---|---:|---|---:|---|---|---|
| pea globulin isolate, lab-made (leg 35 / vic 45 / convic 9 / LOX 7 %) | free SH (1.5 M GdnCl-exposed) | 2.1 ± 0.1 | µmol S−/g protein | **0.0021** | native, 2 wt %, 10 mM phosphate + 5 mM NaCl, pH 7.2 | Table 1a | measured (n = 3) |
| same | S-S | 4.2 ± 0.1 | µmol S−S/g protein | **0.0042** | native | Table 1a | measured |
| same | half-cystine | 2.1 + 2 x 4.2 = 10.5 | not printed | 0.0105 | native | arithmetic | inferred |
| same | free SH after heating | 7.1 ± 0.2 | µmol S−/g protein | **0.0071** | 85 C, 60 min, 2 wt %, pH 7.2, I ~ 0.03 | Table 1b | measured |
| same | S-S after heating | 1.8 ± 0.3 | µmol S−S/g protein | **0.0018** | same | Table 1b | measured |
| BLG, lab-purified from WPI, variant unstated | free SH | 42.5 ± 0.2 | µmol S−/g protein | **0.0425** | native, same medium | Table 1a | measured |
| same | S-S | 102.5 ± 0.3 | µmol S−S/g protein | **0.1025** | native | Table 1a | measured |
| same | half-cystine | 247.5 | not printed | 0.2475 | native | arithmetic | inferred |
| same | free SH after heating | 17.2 ± 0.3 | µmol S−/g protein | **0.0172** | 85 C, 60 min | Table 1b | measured |
| same | S-S after heating | 114.1 ± 0.1 | µmol S−S/g protein | **0.1141** | same | Table 1b | measured |
| βlg/Glob 30/70, 50/50, 70/30 heated | free SH / S-S | 10.7 / 35.8; 12 / 59; 14 / 89 | µmol/g protein | 0.0107 / 0.0358; 0.012 / 0.059; 0.014 / 0.089 | 85 C, 60 min | Table 1b | measured |

**Cross-check against the repo's computed BLG sites** (`protein_matrices.yml`: 1 Cys, 2 S-S, 5
half-cystine per 18 362 Da):

| quantity | computed (mmol/g) | Chihi measured (mmol/g) | measured / computed |
|---|---:|---:|---:|
| free thiol | 0.05446 | 0.0425 | **0.78** |
| disulfide | 0.10892 | 0.1025 | **0.94** |
| half-cystine | 0.2723 | 0.2475 | **0.91** |

The 9 % half-cystine shortfall is of the size expected from the 92 % protein purity being applied on a
Kjeldahl basis that counts non-BLG nitrogen (α-lactalbumin, 8 half-cystine per 14.2 kDa, would push the
other way; peptides and NPN push this way) and from ε = 12 900 vs 13 600 (a 5 % scale factor by
itself). The free-thiol shortfall (78 %) is larger than the S-S shortfall (94 %), consistent with some
Cys121 already oxidised or in non-native dimers in the isolate ("negligible traces of non-native
dimeric βlg" on SDS-PAGE). Verdict: the sequence-computed BLG densities are confirmed to within
10-20 %, with the measured free thiol on the low side.

## 5. Flags

1. **DTNB stoichiometry as printed is insufficient for BLG** (§2): 0.10 µmol DTNB against ~0.55 µmol
   free SH or ~3.2 µmol total SH in the stated 650 µL of 2 wt % protein. Since the BLG results land at
   78-94 % of theory, the samples must have been diluted into the 0-60 µM NAC calibration window; the
   dilution is unstated. The numbers are used as printed; the method description is incomplete.
2. **Pea Glob S-S is low against expectation.** Legumin (35 % of the fraction) is quoted in the paper's
   own introduction as carrying "between two and seven cysteine" per ~60 kDa subunit, i.e. 33-117
   µmol half-cystine per g legumin, or 12-41 µmol/g of this Glob fraction; the measured half-cystine is
   10.5 µmol/g, at or below the bottom of that range. Either the legumin here is Cys-poor, the 20 mM DTT
   / 2 h / 25 C reduction without denaturant was incomplete (GdnCl is added only at the DTNB step, after
   reduction and dialysis), or reoxidation occurred during the "extensive dialysis". The authors call
   the native values "lower and comparable with" Mession et al. 2013. Carry 0.0042 mmol/g as measured
   but flag it as a probable LOWER BOUND on pea-globulin disulfide.
3. **Heated Glob shows S-S falling (4.2 -> 1.8) while free SH rises (2.1 -> 7.1)** with total conserved:
   heating REDUCES net disulfide in pea globulin under these conditions, the opposite of BLG. Directional
   hold-out for any model that oxidises plant-protein thiols on heating. The SDS-PAGE shows ~10 % of
   legumin nonetheless enters new intermolecular S-S, so the Ellman net figure hides an exchange.
4. **BLG variant unstated**; WPI-derived, so mixed A/B. Repo BLG entry is variant A (18 362 Da);
   variant B is 18 276 Da, which changes the computed densities by +0.5 %, immaterial.
5. **"2 wt %" is ambiguous** between powder and protein for the stock solutions; it does not affect
   per-g-protein values.
6. **Unheated mixtures "nd"** for SH / S-S; native mixture values, if needed, are the mass-weighted
   single-protein values (additivity holds for the heated 30/70 and 50/50 within 1 µmol/g; 70/30 S-S
   is 9 µmol/g above additive, flag 4 in §3).
7. **ε = 12 900 M^-1 cm^-1** (NAC-calibrated) vs the usual 13 600 (Ellman) / 14 150 (Riddles): all
   Chihi SH values scale by 0.95 if re-expressed on 13 600.
8. **Ionic strength during heating** is ~0.03 M (10 mM phosphate + 5 mM NaCl), far below food matrices;
   the aggregate-size results (Tables 2, 3) are specific to that.
9. **ASAP version**: no volume/page; cite by DOI.
