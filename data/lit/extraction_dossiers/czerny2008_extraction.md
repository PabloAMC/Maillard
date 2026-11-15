# Czerny et al. 2008 — EXTRACTION (orthonasal odour thresholds of 84 key food odorants in tap water, one Garching (TU Munich) panel of 13-22 assessors, purity-checked references; recognition and detection threshold per compound, in ug/L of water)

### ONE PANEL, ONE METHOD, ONE MEDIUM FOR 84 ODORANTS. The threshold source perception step 7.1 keys first (`data/perception/thresholds.yml`): 15 of the lever table's 48 odorants have a Table 1 row, and a sixteenth (2-ethyl-3,6-dimethylpyrazine) a quoted value in the text.

**Source on disk:** `data/articles/czerny2008.pdf` (9 pp., Eur Food Res Technol 228:265-273). Read in full
from the `pdftotext -layout` text layer (born-digital; Table 1 survives with its columns). The text layer renders
"ug" as `lg`, the multiplication sign as `9` and Greek prefixes as Latin letters; the footnote letter of each
compound (procedure a / b / c) is glued to its name (`Hexanala`). Values read off the layer and checked against
the column positions row by row.

## 0. Identity

| field | value |
|---|---|
| Title | "Re-investigation on odour thresholds of key food aroma compounds and development of an aroma language based on odour qualities of defined aqueous odorant solutions" |
| Authors | Michael Czerny, Martin Christlbauer, Monika Christlbauer, Anja Fischer, Michael Granvogl, Michaela Hammer, Cornelia Hartl, Noelia Moran Hernandez, Peter Schieberle (corresponding) — German Research Centre for Food Chemistry and Department of Chemistry, Technical University of Munich, Garching (as printed) |
| Venue | Eur Food Res Technol **228** (2008) 265-273. Received 20 May 2008, revised 10 July 2008, accepted 17 July 2008, published online 5 August 2008 |
| DOI | 10.1007/s00217-008-0931-x (printed on p. 265) |
| Medium | **tap water** (Table 1 heading: "lg/L of water"); odorants dosed from an ethanolic stock (0.1 mL per L, procedures A and C) or directly (procedure B) |
| Route | **orthonasal** (section heading "Determination of orthonasal odour thresholds in water"; vessels opened and sniffed) |
| Quantity | **recognition** and **detection** odour threshold per compound; panel value = geometric mean of individual thresholds |

## 1. Method (what a criterion field must say)

- **Purity check first.** Every reference was run by HRGC-O / aroma extract dilution analysis; used only if the
  target's FD factor was at least 100x that of the most intense odour-active contaminant (p. 268). (E)-2-decenal
  failed and was purified by preparative GC; octanal was found degraded to 28 % after 15 months at -30 C and
  re-ordered (p. 270).
- **Series.** Stock about 50x the estimated recognition threshold, diluted 1:3 in tap water. Triangular series
  (one odorant sample, two tap-water blanks, 20 mL in 45 mL Teflon vessels), presented in **decreasing**
  concentration with the sample position rotated between assessors (p. 269).
- **Individual threshold** OT_i = sqrt(C_x C_(x-1)) in the definitions (C_x the lowest concentration correctly
  selected, C_(x-1) the highest incorrectly selected; the displayed formula prints the index as x+1, a typo
  against its own definitions). **Panel threshold** = geometric mean of the OT_i.
- **Panel.** 13-22 trained assessors present on the day (of up to 22, aged 25-39), weekly training; 21 +/- 1 C
  sensory room. No chance correction is stated; the triangle design makes the guess rate 1/3.
- **Recognition vs detection.** Detection = a difference perceived; recognition = the odour quality named. The
  ratio recognition / detection lies between 1.3 (dimethyl sulphide) and 7.3 (1-butanol) (p. 270).

**Keyed criterion:** the **detection** threshold, as `detection_triangle_series_geometric_mean`. One criterion
per odorant; the recognition value of the same row is quoted in each record's notes, not keyed.

## 2. Table 1 rows used (lever-table odorants), as printed

`No. | compound (procedure) | odour quality | recognition | detection | literature range | refs`, ug/L of water.

| No. | compound as printed | quality | recognition | detection | literature range |
|---:|---|---|---:|---:|---|
| 1 | Acetaldehyde (a) | Fresh, green | 63 | 25 | 15-120 |
| 8 | 2,3-Butanedione (a) | Buttery | 2.9 | 1.0 | 3-6.5 |
| 13 | (E,E)-2,4-Decadienal (a) | Fatty, deep-fried | 0.077 | 0.027 | 0.07-0.2 |
| 19 | Dimethyl trisulphide (a) | Cabbage-like | 0.016 | 0.0099 | 0.01 |
| 34 | 2-Furfurylthiol (a) | Sulphury, burnt | 0.080 | 0.036 | 0.006-0.047 |
| 37 | Hexanal (a) | Green, grassy | 10 | 2.4 | 4.5-479 |
| 43 | 4-Hydroxy-2,5-dimethyl-3(2H)-furanone (a) | Caramel-like | 160 | 40 | 3.2-1,700 |
| 50 | (R/S)-2-Methylbutanal (a) | Malty | 4.4 | 1.5 | 3-13 |
| 51 | 3-Methylbutanal (a) | Malty | 1.2 | 0.50 | 0.2-170 |
| 56 | Bis-(2-methyl-3-furyl) disulphide (a) | Meaty | 0.00076 | 0.00032 | 0.00002-0.00008 |
| 61 | Methylpropanal (a) | Malty | 1.9 | 0.49 | 1-10 |
| 64 | 3-(Methylthio)-propanal (a) | Cooked potato-like | 1.4 | 0.43 | 0.2-1.8 |
| 70 | Nonanal (a) | Citrus-like, soapy | 8.0 | 2.8 | 1-98 |
| 71 | (E)-2-Nonenal (a) | Fatty, green | 0.69 | 0.19 | 0.08-0.4 |
| 74 | Octanal (a) | Citrus-like, green | 6.9 | 3.4 | 0.58-6.4 |

All fifteen are procedure A (polar; ethanolic stock into tap water). Name mappings to `data/keys/compounds.yml`:
row 19 "trisulphide" = dimethyl trisulfide; row 43 = DMHF / Furaneol (`hdmf`); row 50 the racemate of
2-methylbutanal; row 61 = 2-methylpropanal; row 64 = methional.

**In the text, not Table 1 (p. 270):** Guadagni et al. [22] measured 0.8 ug/L of water for
2-ethyl-3,6-dimethylpyrazine; finding traces of the more potent 2-ethyl-3,5-dimethylpyrazine in their sample,
they later reported **8.6 ug/L** after purification. Keyed as a secondary quotation (Czerny give no reference for
the re-examination).

**Lever-table odorants with no Czerny value:** MFT, phenylacetaldehyde, pentanal, heptanal, (E)-2-octenal,
1-octen-3-ol (1-octen-3-one only, row 75), 1-hexanol, HEMF, norfuraneol, hydrogen sulphide, every alkylpyrazine
of the table, the furans, furfurals and furfuryl disulphides, benzaldehyde, 2-acetylthiazole, 2-pentylpyridine.

## 3. Bis(2-methyl-3-furyl) disulfide: the 30 000x discrepancy (SCOPING_2026-11.md sec. 2)

| source | value as printed | in ug/L | cited? |
|---|---|---:|---|
| Czerny 2008 Table 1 row 56, detection | 0.00032 ug/L | 3.2e-4 | measured here |
| Czerny 2008 Table 1 row 56, recognition | 0.00076 ug/L | 7.6e-4 | measured here |
| Czerny 2008 Table 1 row 56, literature range | 0.00002-0.00008 | 2e-5 - 8e-5 | refs [61, 62] |
| Kerler 2010 Table 3.3 row 10 | 0.00002 ug/kg | 2e-5 | Buttery et al. 1984 |
| Thong 2024 Table 4 row 8 | 0.0007 ppt | 7e-7 | no |
| `desirable_targets.yml` | 0.02 ug/kg | 2e-2 | no |

Thong's "0.0007 ppt" matches Czerny's recognition threshold 0.00076 ug/L (0.0007 ppb, truncated) to the digit if
**ppb was misprinted as ppt**; read as ppb, it is 0.7 ng/L, inside the published spread. That reading is an
inference (Thong cite nothing), but the coincidence of digits with the 2008 panel's value and the 1000x size of the
slip make it the likely explanation (my estimate: about 80 %). On it, the 30 000x gap splits into a 1000x unit
slip and a real but ordinary spread: `desirable_targets.yml`'s 0.02 ug/kg is 63x Czerny's detection threshold
(26x the recognition one), and the measured values themselves span 2e-5 to 7.6e-4 ug/L, 1.6 dex, between
laboratories. The uncited file's number is still the outlier (1.4 - 3 dex above every measured value) and is not
used.

## 4. Use in the model

- **Keyed:** the 15 detection thresholds above plus the 2-ethyl-3,6-dimethylpyrazine quotation, into
  `data/perception/thresholds.yml` (medium water, ug/L, orthonasal, 21 C). Reference inputs only; nothing in
  `src/kinetic_core/` reads them (the engine is paused).
- **Not keyed:** recognition thresholds (one criterion per odorant); the "literature range" column (other
  laboratories, unnamed methods); the 69 Table 1 compounds outside the lever table.
- **What it does not fix:** water is not a plant-protein product. Matrix suppression in a fat- or protein-rich
  product is 0 to +2 dex on the effective threshold (SCOPING sec. 4, row 7.3), so a water threshold never makes a
  "stays below" call in such a product.

## 5. Flags

1. The OT_i formula's index typo (x+1 printed, x-1 defined); the definitions are taken as meant.
2. Ethanol carrier at 0.01 % v/v in procedures A and C; its own odour is far below threshold (ethanol row 25:
   detection 990 000 ug/L), so no correction is needed.
3. Thresholds are panel geometric means with no dispersion printed; the between-session repeat of Table 2 (three
   acids only) differs by up to 5x (phenylacetic acid detection 6,100 vs 1,200 ug/L), a rough size for the
   within-laboratory spread.
