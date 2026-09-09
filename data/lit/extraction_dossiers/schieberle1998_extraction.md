# Schieberle & Hofmann 1998 — EXTRACTION (dry-heated 180 °C / 6 min versus aqueous 145 °C / 20 min cysteine + ribose and cysteine + rhamnose, SIDA levels of the key odorants)
### The dry-versus-aqueous companion of Hofmann 1998 (same laboratory, same pots); its aqueous column is Hofmann 1998's ribose and rhamnose rows reprinted, its new content is the mercaptoketone, thiazoline, thiazine, furfural, furanone and pyrazine levels in both regimes.

**Source on disk:** `data/articles/schieberle1998.pdf` (owner's download, 2026-09-08; 11 pages, scanned ACS
chapter with an OCR text layer). Read from the pdftotext extraction in the scratchpad; the three
quantitative tables (IV, VI, VII) were re-checked against page renders (pdftoppm, 110 dpi) because the
OCR garbled their unit headers: Table IV and VI read "µg/100 mmol ribose (rhamnose)" and Table VII reads
"Amount (µg)" on the page (OCR gave "lig", "jug" and "ng"). Figures 1 to 5 (labelled standards, purge-and-trap
scheme, two spider webs, one mechanism) carry no numbers the text does not also give and are FIGURE-ONLY.

## 0. Identity

| field | value |
|---|---|
| Title | "Characterization of Key Odorants in Dry-Heated Cysteine–Carbohydrate Mixtures: Comparison with Aqueous Reaction Systems" |
| Authors | P. Schieberle and T. Hofmann, Deutsche Forschungsanstalt für Lebensmittelchemie, Garching |
| Venue | ACS Symposium Series 705, *Flavor Analysis* (Mussinan & Morello, eds.), 1998, **Chapter 26**, pp. 320-330 |
| DOI | 10.1021/bk-1998-0705.ch026 (printed on every page; publication date 10 Sept 1998). The ACS download file name carries "ch032"; the chapter number printed in the PDF is 26 — cite the DOI, not the file name. |
| Sibling papers | ref 3 = Hofmann & Schieberle 1995 JAFC 43, 2187 (aqueous ribose/cysteine AEDA; on disk `hofmann1995.pdf`); ref 4 = Hofmann & Schieberle 1997 JAFC 45, 898 (glucose, rhamnose AEDA; not on disk); "Hofmann and Schieberle, J. Agric. Food Chem., 1998, submitted" = the full paper behind Tables II, V and the labelled standards (probably `hofmann1998.pdf`, JAFC 46, 235, whose Table 2 dry column is the same run — see §4.4); ref 6 = Schieberle & Hofmann 1996 (Flavour Science, RSC) on MFT losses during work-up |
| Repo cross-reference | `hofmann1998_reconciliation.md` (the JAFC paper: 145 °C aqueous / 180 °C dry, FFT and MFT only); `schieberle2000_extraction.md` (Table III there reprints this chapter's aqueous column with two more odorants; Table IV there is the 100 °C time course) |

## 1. Why it matters

The sulfur lane's formation and sink constants are anchored on Hofmann 1998's 145 °C / 20 min aqueous pot,
and that paper prints only FFT and MFT. This chapter is the same laboratory, the same year and the same
two pots (cysteine 3.33 mmol + ribose 10 mmol; the aqueous one in 100 mL pH-5 phosphate, the dry one on
silica gel), but quantifies eight further odorants by stable isotope dilution in both regimes: the
mercaptoketone 3-mercapto-2-pentanone, 2-acetyl- and 2-propionyl-2-thiazoline, 5-acetyl-2,3-dihydro-1,4-
thiazine, furan-2-aldehyde and three pyrazines (ribose), and 5-methyl-2-furfurylthiol, HDMF, the pyranone
and 5-methylfurfural (rhamnose). For the sink re-fit the useful content is (i) the 145 °C aqueous levels of
3-mercapto-2-pentanone and the thiazine, which are the 145 °C anchors of the 100 °C time series in
Schieberle 2000 Table IV and are printed nowhere in Hofmann 1998, and (ii) the within-study dry/aqueous
ratios, where the thiols rise (FFT ×8, MFT ×1.3) while the mercaptoketone falls ×6 and the thiazine
falls ×42. The regime jump is confounded (water, temperature, time, buffer strength all change at once,
see `hofmann1998_reconciliation.md` §1.2), so these are direction-and-magnitude hold-out shapes, not an
a_w axis.

## 2. Methods as they matter to a model

- **Dry system ("I").** Table IV footnote c: "Cysteine (3.33 mmol) and ribose (10 mmol) were intimately
  mixed with silica gel (2.7 g containing 0.3 g of the phosphate buffer) and were reacted for 6 min at
  180 °C." Text: "mixed with silica gel (3 g) containing 10 % of an aqueous sodium phosphate solution
  (0,5 mol/L; pH 5.0) and heated for 6 min at 180 °C in a closed vessel. Phosphate was used as a catalyst
  but not to stabilize the pH." So: ~0.3 g of liquid on ~3 g of solid; the charge is 10 mmol ribose
  (1.50 g) + 3.33 mmol cysteine (0.40 g) + 3 g silica/buffer, i.e. ~0.3 g water in ~4.9 g total, ~6 % w/w.
  No water activity, no vessel type, no heating device (oven vs block) and no headspace volume are
  stated. ⚠ Buffer molarity conflicts between sources for the same arm: 0.5 mol/L here (text), "0.3 g of
  the phosphate buffer" (footnote, molarity by reference to footnote d = 0.5 mol/L), but 300 µL of
  **2 mol/L** in Hofmann 1998 JAFC (Experimental, p. 236). The water amount agrees (0.3 g ≈ 300 µL); the
  phosphate load does not (0.15 vs 0.6 mmol). Record as unresolved.
- **Aqueous system ("II").** Table IV footnote d: "Cysteine (3.33 mmol) and ribose (10 mmol) were dissolved
  in phosphate buffer (100 mL; 0.5 mol/L; pH 5.0) and reacted for 20 min at 145 °C." ⚠ The running text
  says "100 mL of phosphate buffer; 0.1 mol/L; pH 5.0" for the same previous studies; Hofmann 1998 JAFC
  and Schieberle 2000 both say 0.5 mol/L; read 0.1 as a slip. Vessel not named here (JAFC: 200 mL Roth
  Type II laboratory autoclave). Concentrations 100 mmol/L ribose, 33.3 mmol/L cysteine, cys:sugar 1:3.
- **Rhamnose pots**: same charges and conditions with rhamnose in place of ribose (Table VI, "for
  footnotes see Table IV").
- **Ten-fold batch (Table VII)**: "Cysteine (33 mmol) and the carbohydrate (100 mmol) were mixed with
  silica gel (27 g containing 3 g of phosphate buffer [0.5 mol/L; pH 5.0]) and were reacted for 6 min
  at 180 °C." Ribose, rhamnose and glucose.
- **Isolation for AEDA**: diethyl ether extraction, then "sublimation in vacuo" (high-vacuum
  distillation, ref 3); HRGC/olfactometry with stepwise dilution; identification against reference
  compounds.
- **Quantification**: stable isotope dilution assays with deuterium- or carbon-13-labelled internal
  standards (structures in Figure 1, FIGURE-ONLY; syntheses in the "submitted" JAFC paper). Because
  MFT and its labelled standard "are completely degraded" when an extraction/distillation/concentration
  work-up is used (ref 6), the thiols were enriched by **purge and trap** (TCT system, Chrompack):
  helium flushing of the sample in 15 mL of 0.5 mol/L phosphate pH 5.0 for 20 min onto a trap, then
  GC-MS with selected-ion mass chromatography. Recoveries (Table III, 2-5 µg of each thiol spiked)
  95-104 %, "differed not more than 5 % from the actual amounts". Which analytes beyond the seven
  thiols of Table III had labelled standards is not stated in this chapter. Replicates: not stated
  here (Hofmann 1998 JAFC: triplicates, ±10 %).
- **Reporting unit**: "The concentrations displayed in Table IV are based on 10 mmol of ribose; the odor
  activity values are based on dissolving each of the processed mixtures in 100 mL of water." The table
  header reads **µg/100 mmol ribose**. The two are reconciled by the OAV arithmetic: a 10-mmol pot
  dissolved in 100 mL gives µg/L = 10 × (µg per pot) = µg per 100 mmol ribose, and every printed OAV is
  exactly conc/threshold (972/0.01 = 97200). So **µg per 100 mmol carbohydrate = µg/L of the 100 mL pot;
  divide by 10 for µg per the standard pot (10 mmol sugar, 3.33 mmol cysteine)**. Hofmann 1998 JAFC
  prints the per-pot values (FFT 12.1, 97.2; MFT 19.8, 25.1 µg) and they match this chapter's 121, 972,
  198, 251 exactly (§4.4).
- **Odor thresholds**: in water, triangle test, seven-member panel (Table IV footnote b).
- **Sensory profiles**: ten-member panel, spider webs (Figures 3, 4; FIGURE-ONLY). Text: both ribose
  systems "roasty and meat-like" predominated; "an earthy note was only detectable in the dry-heated
  mixture". Rhamnose: caramel and seasoning notes in the aqueous mixture, roasty (5-methyl-2-
  furfurylthiol) in the dry one.

## 3. Tables re-typed

### Table I. "Most intense odorants (FD ≥ 256 in at least one mixture) generated by thermal treatment of cysteine (C) and carbohydrates in aqueous solution (ribose: Rib, rhamnose: Rha, glucose: Glc). Data from (3,4)"

| odorant | odor quality | FD C/Rib | FD C/Rha | FD C/Glc |
|---|---|---:|---:|---:|
| 2-furfurylthiol | roasty, coffee-like | 1024 | 512 | 1024 |
| 3-mercapto-2-pentanone | catty, sulfury | 512 | 128 | 512 |
| 2-methyl-3-furanthiol | meat-like | 256 | <4 | >4 (as printed; probably <4) |
| 5-acetyl-2,3-dihydro-1,4-thiazine | roasty, popcorn | 256 | 512 | 1024 |
| 3-mercaptobutanone | sulfury | 128 | 32 | 512 |
| 2-(1-mercaptoethyl)furan | burnt | <1 | <1 | 256 |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone | caramel-like | 32 | 65536 | 512 |
| 5-methyl-2-furfurylthiol | roasty, coffee-like | <4 | 2048 | <4 |
| 3-hydroxy-6-methylpyran-2-one | seasoning-like | <4 | 16384 | <4 |

### Table II. "Most odor-active volatiles identified in the dry-heated cysteine/ribose mixture" (11 of the 24 odorants in the FD range 2-16384)

| odorant | FD | odorant | FD |
|---|---:|---|---:|
| 2-furfurylthiol | 16384 | 2,3-diethyl-5-methylpyrazine | 128 |
| 2-acetyl-2-thiazoline | 1024 | 2-ethenyl-3,5-dimethylpyrazine | 128 |
| 2-methyl-3-furanthiol | 256 | 4-hydroxy-2,5-dimethyl-3(2H)-furanone | 128 |
| furan-2-aldehyde | 256 | 3-mercapto-2-pentanone | 64 |
| 2-ethyl-3,5-dimethylpyrazine | 256 | 3-hydroxy-4,5-dimethyl-2(5H)-furanone (sotolon) | 64 |
| 2-propionyl-2-thiazoline | 256 | | |

Compared with Table I (aqueous, FFT 1024, MFT 256, MP 512, thiazine 256): FFT ×16 in FD, MFT unchanged,
3-mercapto-2-pentanone ÷8, thiazine off the list. Authors: "especially the FD factors of 2-furfurylthiol and
2-acetyl-2-thiazoline as well as those of the three pyrazines were increased, whereas the FD-factors of the
popcorn-like 5-acetyl-2,3-dihydro-1,4-thiazine and the sulfury 3-mercapto-2-pentanone were decreased in the
dry-heated system."

### Table III. "Recoveries of 7 thiols determined by stable isotope dilution assays and using the 'purge and trap' technique as the enrichment procedure"

| thiol | recovery (%) |
|---|---:|
| 3-mercapto-2-butanone | 99 |
| 3-mercapto-2-pentanone | 95 |
| 2-furfurylthiol | 96 |
| 2-(1-mercaptoethyl)furan | 97 |
| 2-methyl-3-furanthiol | 104 |
| 2-thiophenemethanethiol | 96 |
| 2-(1-mercaptoethyl)thiophene | 98 (OCR "9S") |

Footnote a: 2-5 µg of each thiol and 2-5 µg of each labelled standard in 15 mL 0.5 mol/L phosphate pH 5.0,
swept onto the trap with helium for 20 min. Footnote b: SIDA by mass chromatography of selected ions.

### Table IV. "Concentrations and odor activity values (OAV) of selected key odorants generated from ribose and cysteine — Influence of the reaction conditions"

Units: conc. in **µg/100 mmol ribose** (= µg/L of the 100 mL pot; ÷10 = µg per pot); OAV = conc/threshold.
I = dry-heated (footnote c: 180 °C, 6 min, silica gel); II = aqueous (footnote d: 145 °C, 20 min, 100 mL
0.5 mol/L phosphate pH 5.0). Odor threshold in µg/L water in parentheses (footnote b).

| odorant (threshold µg/L) | conc. I dry | conc. II aqueous | OAV I | OAV II |
|---|---:|---:|---:|---:|
| 2-furfurylthiol (0.01) | 972 | 121 | 97200 | 12100 |
| 2-acetyl-2-thiazoline (1.0) | 49 | 7 | 49 | 7 |
| 2-propionyl-2-thiazoline (1.0) | 18 | <1 | 18 | <1 |
| 2-ethyl-3,5-dimethylpyrazine (0.1) | 11 | <0.1 | 110 | <1 |
| 2,3-diethyl-5-methylpyrazine (0.1) | 3 | <0.1 | 30 | <1 |
| 2-ethenyl-3,5-dimethylpyrazine (0.1) | 5 | <0.1 | 50 | <1 |
| furan-2-aldehyde (12000) | 79000 | 52 | 7 | <1 |
| 5-acetyl-2,3-dihydro-1,4-thiazine (1.25) | 10 | 424 | 8 | 340 |
| 3-mercapto-2-pentanone (0.7) | 101 | 599 | 144 | 856 |
| 2-methyl-3-furanthiol (0.007) | 251 | 198 | 35857 | 28286 |

Authors: "The quantitative data corroborated the results of the AEDA and indicated an increase of
2-furfurylthiol, 2-acetyl-2-thiazoline and three pyrazines in the dry-heated (I) compared with the aqueous
system (II). On the other hand, the amounts of 5-acetyl-2,3-dihydro-1,4-thiazine and 3-mercapto-2-pentanone
were significantly decreased." And: "in the dry-heated mixture FFT is higher than MFT, whereas the reverse
was true for the aqueous system." On furfural: "Although present in high concentrations, due to its high
odor threshold, furan-2-aldehyde showed a very low odor contribution to the dry-heated system."

### Table V. "Most odor-active volatiles identified in the dry-heated cysteine/rhamnose mixture" (12 of 22 odorants, FD range 2-4096)

| odorant | FD | odorant | FD |
|---|---:|---|---:|
| 5-methyl-2-furfurylthiol | 4096 | 5-methylfuran-2-aldehyde | 512 |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone (Furaneol) | 2048 | 2-ethyl-3,5-dimethylpyrazine | 256 |
| (Z)-2-propenyl-3,5-dimethylpyrazine | 1024 | 2,3-diethyl-5-methylpyrazine | 128 |
| 2-acetyl-2-thiazoline | 1024 | unknown (sweet, roasty) | 128 |
| 2-propionyl-2-thiazoline | 1024 | 2-furfurylthiol | 32 |
| | | 3-hydroxy-6-methyl-2(2H)-pyranone | 32 |
| | | 3-hydroxy-4,5-dimethyl-2(5H)-furanone | 32 |

Compared with Table I (aqueous rhamnose: HDMF 65536, pyranone 16384, 5-MFFT 2048): "a drastic decrease in the
FD-factors of 4-hydroxy-2,5-dimethyl-3(2H)-furanone (caramel-like) and 3-hydroxy-6-methyl-2(2H)-pyranone
(seasoning-like) was detected."

### Table VI. "Concentrations and odor activity values (OAVs) of selected key odorants generated from rhamnose and cysteine — Influence of the reaction conditions" (units and footnotes as Table IV, rhamnose in place of ribose)

| odorant (threshold µg/L) | conc. I dry | conc. II aqueous | OAV I | OAV II |
|---|---:|---:|---:|---:|
| 5-methyl-2-furfurylthiol (0.05) | 1156 | 65 | 24083 | 1354 |
| 2-acetyl-2-thiazoline (1.0) | 48 | 6 | 48 | 6 |
| 2-propionyl-2-thiazoline (1.0) | 18 | <1 | 18 | <1 |
| 2-ethyl-3,5-dimethylpyrazine (0.1) | 12 | <0.1 | 120 | <1 |
| (Z)-2-propenyl-3,5-dimethylpyrazine (0.1) | 3 | <0.1 | 30 | <1 |
| 2,3-diethyl-5-methylpyrazine (0.1) | 1 | <0.1 | 10 | <1 |
| 5-methylfuran-2-aldehyde (4500) | 147000 | <10 | 331 | <1 |
| 4-hydroxy-2,5-dimethyl-3(2H)-furanone (10) | 35600 | 198000 | 3560 | 19800 |
| 3-hydroxy-6-methyl-2(2H)-pyranone (15) | 1355 | 245300 | 90 | 16353 |
| 2-methyl-3-furanthiol (0.007) | 31 | 8 | 4429 | 1143 |
| 2-furfurylthiol (0.01) | 4 | 8 | 400 | 800 |

Note the printed OAV for 5-MFFT dry (24083 ≠ 1156/0.05 = 23120) and for 5-methylfurfural (331 ≠ 147000/4500 =
32.7) do not reproduce from the printed concentration and threshold; all other cells do. Use the
concentrations, not those two OAVs.

### Table VII. "Amounts of odor-active pyrazines formed in cysteine/carbohydrate mixtures" (ten-fold dry batch: cysteine 33 mmol + carbohydrate 100 mmol on 27 g silica gel + 3 g of 0.5 mol/L phosphate pH 5.0; 6 min, 180 °C). Unit on the page: **µg** (OCR "ng").

| pyrazine | ribose | rhamnose | glucose |
|---|---:|---:|---:|
| 2-ethyl-3,5-dimethyl- | 11 | 12 | 15 |
| 2-ethenyl-3,5-dimethyl- | 5 | <0.1 | <0.1 |
| (Z)-2-propenyl-3,5-dimethyl- | <0.1 | 3 | <0.1 |
| 2,3-diethyl-5-methyl- | 3 | 1 | 2 |

The ribose and rhamnose columns equal Table IV/VI column I cell for cell, which means the pyrazine values
in Tables IV and VI are these 100-mmol-batch amounts (µg per 100 mmol sugar) rather than a separate
10-mmol run; consistent with the unit reading in §2. Authors: EDMP "is formed from the three sugars in
nearly equal amounts"; the ethenyl- and propenylpyrazine "were formed exclusively from either ribose or
rhamnose, respectively"; the four pyrazines "are predominantly formed in the dry-heated systems".

## 4. Numbers the repository can use

Conversion: with X = µg per 100 mmol carbohydrate, µg per standard pot = X/10; mol % of carbohydrate =
(X/MW)/100 000 µmol × 100 = X/(MW × 1000); mol % of cysteine (33.3 mmol per 100 mmol sugar) = 3.0 × that.
Molar masses used: FFT and MFT C5H6OS 114.17; 3-mercapto-2-pentanone C5H10OS 118.20; 2-acetyl-2-thiazoline
C5H7NOS 129.18; 2-propionyl-2-thiazoline and 5-acetyl-2,3-dihydro-1,4-thiazine C6H9NOS 143.21; furfural
C5H4O2 96.08; 2-ethyl-3,5-dimethylpyrazine C8H12N2 136.19; 5-methyl-2-furfurylthiol C6H8OS 128.19;
5-methylfurfural C6H6O2 110.11; HDMF C6H8O3 128.13; 3-hydroxy-6-methyl-2H-pyran-2-one C6H6O3 126.11.
Worked example: FFT dry 972 µg / 114.17 = 8.51 µmol per 100 mmol ribose = 0.00851 mol % of ribose
(0.0256 mol % of cysteine); aqueous 121/114.17 = 1.06 µmol = 0.00106 mol %.

### 4.1 Levels (validate) and within-study ratios (may be fitted), ribose pots

| quantity | dry 180 °C / 6 min | aqueous 145 °C / 20 min | dry/aq ratio | evidence class | source |
|---|---:|---:|---:|---|---|
| FFT, µg per pot / mol % ribose | 97.2 / 0.00851 | 12.1 / 0.00106 | 8.03 | measured_level (already held from Hofmann 1998 T2) | Table IV |
| MFT, µg per pot / mol % ribose | 25.1 / 0.00220 | 19.8 / 0.00173 | 1.27 | measured_level (already held) | Table IV |
| 3-mercapto-2-pentanone, µg per pot / mol % | 10.1 / 0.000854 | **59.9 / 0.00507** | **0.169** | measured_level, **new 145 °C anchor** for the Schieberle 2000 100 °C series (2.1/10.5/79/85 µg) | Table IV |
| 5-acetyl-2,3-dihydro-1,4-thiazine, µg per pot / mol % | 1.0 / 0.0000698 | **42.4 / 0.00296** | **0.024** | measured_level, new 145 °C anchor (Schieberle 2000: peaks 5.2 µg at 60 min / 100 °C, <0.1 at 12 h) | Table IV |
| 2-acetyl-2-thiazoline, µg per pot / mol % | 4.9 / 0.000379 | 0.7 / 0.0000542 | 7.0 | measured_level | Table IV |
| 2-propionyl-2-thiazoline, µg per pot | 1.8 | <0.1 | >18 | measured_level | Table IV |
| furan-2-aldehyde, µg per pot / mol % | 7900 / **0.822** | 5.2 / 0.000541 | ~1500 | measured_level; aqueous cell conflicts with Hofmann 1998 T5 (§5.2) | Table IV |
| 2-ethyl-3,5-dimethylpyrazine, µg per pot / mol % | 1.1 / 0.0000808 | <0.01 | >110 | measured_level (from the 100-mmol batch) | Tables IV, VII |
| MFT/FFT mass ratio | 0.258 | 1.64 | — | within_study_ratio | Table IV |
| MP/MFT mass ratio | 0.40 | 3.03 | — | within_study_ratio | Table IV |
| thiazine/MFT mass ratio | 0.040 | 2.14 | — | within_study_ratio | Table IV |
| 3-mercapto-2-pentanone at 100 °C 6 h ÷ 145 °C 20 min (Schieberle 2000 Table IV 79 µg ÷ 59.9) | — | — | 1.32 | within_study_ratio across the two chapters (same pot, same lab) | this Table IV + Schieberle 2000 Table IV |

### 4.2 Levels and ratios, rhamnose pots

| quantity | dry | aqueous | dry/aq | evidence class | source |
|---|---:|---:|---:|---|---|
| 5-methyl-2-furfurylthiol, µg per pot / mol % rhamnose | 115.6 / 0.00902 | 6.5 / 0.000507 | 17.8 | measured_level | Table VI |
| MFT, µg per pot / mol % | 3.1 / 0.000272 | 0.8 / 0.0000701 | 3.9 | measured_level (already held from Hofmann 1998 T2) | Table VI |
| FFT, µg per pot / mol % | 0.4 / 0.000035 | 0.8 / 0.0000701 | 0.5 | measured_level (already held) | Table VI |
| 2-acetyl-2-thiazoline, µg per pot | 4.8 | 0.6 | 8 | measured_level | Table VI |
| HDMF, µg per pot / mol % | 3560 / 0.278 | 19800 / **1.55** | 0.18 | measured_level | Table VI |
| 3-hydroxy-6-methyl-2H-pyran-2-one, µg per pot / mol % | 135.5 / 0.0107 | 24530 / **1.95** | 0.0055 | measured_level | Table VI |
| 5-methylfurfural, µg per pot / mol % | 14700 / 1.34 | <1 | >14000 | measured_level | Table VI |

### 4.3 Directional statements worth holding (short quotes)

- Dry heating moves sulfur from the mercaptoketone and the thiazine into the furanthiols and thiazolines:
  MP ÷6, thiazine ÷42, FFT ×8, MFT ×1.3, 2-AT ×7. A model in which 3-mercapto-2-pentanone is a
  sink-side product of the same H2S pool that feeds MFT and FFT should reproduce the sign flip.
- "5-acetyl-2,3-dihydro-1,4-thiazine ... drastically decreased" under dry heating (Schieberle 2000 repeats
  it); together with its 100 °C peak-and-vanish (Schieberle 2000 Table IV) the thiazine is the one
  odorant in this pot that behaves like a labile intermediate at every condition.
- Furfural is 0.8 mol % of ribose in the dry pot and 0.0005 mol % in the aqueous pot: the dry regime is a
  furfural-rich regime, which is where the FFT route (furfural + H2S, Hofmann 1998 Table 3) would be
  expected to win; the observed FFT ×8 is consistent with that and MFT ×1.3 with the NF route being
  water-dependent.
- Rhamnose: HDMF and the pyranone fall 5× and 180× in the dry pot while 5-methylfurfural appears at
  1.3 mol %; 5-MFFT rises 18×.

### 4.4 Same experiment or new run? (against Hofmann 1998 JAFC and Schieberle 2000)

| odorant | this chapter, aqueous (µg/100 mmol) | Hofmann 1998 JAFC T1/T2 (µg per pot) ×10 | this chapter, dry | JAFC T2 dry ×10 |
|---|---:|---:|---:|---:|
| FFT, ribose | 121 | 121 | 972 | 972 |
| MFT, ribose | 198 | 198 | 251 | 251 |
| FFT, rhamnose | 8 | 8 | 4 | 4 |
| MFT, rhamnose | 8 | 8 | 31 | 31 |

All eight cells agree to the last digit: the FFT and MFT numbers are **the same measurements reprinted at
ten-fold scale**, not a replicate run; do not enter them a second time. Schieberle 2000 Table III (145 °C,
20 min, "based on 33 mmol of cysteine and 100 mmol of the carbohydrate") reprints this chapter's aqueous
column again (FFT 121/8, MP 599/73, MFT 198/8, HDMF and pyranone as here, 5-MFFT <0.1/65, thiazine
424/401) and adds 3-mercapto-2-butanone 342 (ribose) and 141 (rhamnose) µg/100 mmol. ⚠ Note for the
inventory: `hofmann1998_reconciliation.md` records "342 / 200 ppb" as values REFUSED because they appear
nowhere in Hofmann 1998 JAFC; the 342 does exist in this family of papers, as **3-mercapto-2-butanone,
µg/100 mmol ribose (= 34.2 µg per pot = 342 µg/L), Schieberle 2000 Table III**, attributed to the wrong
paper and probably the wrong analyte. The 200 is not found in either chapter (MFT is 198).

## 5. Flags

1. **Unit header.** OCR renders it "lig/100 mmol"; the page reads µg/100 mmol carbohydrate, and the OAV
   arithmetic and the JAFC cross-check both confirm value/10 = µg per 10-mmol pot. Any ingestion must
   apply the ÷10 (or treat the values as µg/L of the 100 mL pot).
2. **Furfural, aqueous pot: 52 µg/100 mmol here (5.2 µg per pot) versus 67.5 µg per pot in Hofmann 1998
   Table 5 ("FA")**, a 13-fold disagreement for what should be the same pot, same lab, same year. No
   explanation in either paper (different standard, different work-up, or a typo). Do not use the aqueous
   furfural level from either source without a third measurement; the dry-pot 79000 has no counterpart.
3. **Buffer strength of the dry arm** is 0.5 mol/L here and 2 mol/L in the JAFC paper; the aqueous buffer is
   0.5 mol/L in the footnote and 0.1 mol/L in the running text. Water amount (0.3 g on 3 g silica) is
   consistent across all descriptions.
4. **The regime jump is four-fold confounded** (a_w, 145→180 °C, 20→6 min, buffer). Ratios in §4.1 are
   hold-out shapes for a dry-regime run of the model, not coefficients for a water term (the existing
   FIT_HOLDOUT_DECLARATION line for Hofmann 1998 T2 dry rows already says HOLD-OUT; the new MP, thiazine
   and 2-AT ratios belong under the same entry).
5. **No replicate statement, no error bars** in this chapter; inherit the JAFC "triplicates, ≤10 %" only
   for the FFT/MFT cells that are demonstrably the same data. Recovery of the purge-and-trap SIDA is
   95-104 % for seven thiols (Table III), the strongest method statement in the sulfur corpus.
6. **Two printed OAVs in Table VI do not reproduce** (5-MFFT dry, 5-methylfurfural dry); the
   concentrations are internally consistent with the JAFC paper, so keep the concentrations.
7. **Table I anomaly**: MFT in the glucose column is printed ">4" where the pattern demands "<4".
8. **Registry keys** (`data/keys/compounds.yml`): MFT `2_methyl_3_furanthiol`; FFT `2_furfurylthiol`;
   bis(2-methyl-3-furyl) disulfide `bis_2_methyl_3_furyl_disulfide`; furfural `furfural`; HDMF `hdmf`;
   2-ethyl-3,5-dimethylpyrazine `2_ethyl_3_5_dimethylpyrazine`; norfuraneol `norfuraneol`. **No key exists**
   for 3-mercapto-2-pentanone, 3-mercapto-2-butanone (= 2-mercapto-3-butanone), 2-acetyl-2-thiazoline,
   2-propionyl-2-thiazoline, 5-acetyl-2,3-dihydro-1,4-thiazine, 5-methyl-2-furfurylthiol,
   3-hydroxy-6-methyl-2H-pyran-2-one, 5-methylfurfural, sotolon, 2-furfuryl methyl disulfide, or the
   ethenyl/propenyl/diethyl pyrazines. The 145 °C anchors in §4.1 for MP and the thiazine therefore need
   keys before they can be ingested.
9. **Sotolon** (3-hydroxy-4,5-dimethyl-2(5H)-furanone) is detected in both dry pots (FD 64 ribose, 32
   rhamnose) but not quantified.
