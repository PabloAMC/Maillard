# Li, Dai, Mao, An, Bai & Kaur 2026 — EXTRACTION (enzymatic pretreatment of wheat gluten and the volatiles of a high-moisture SPI:gluten extrudate)

**Source on disk:** `data/articles/li2026.pdf` (2.0 MB; downloaded 2026-09-11 at the 2026-09-11
reading-list row's request). Read 2026-09-11 via `pdftotext -layout`; Table 2 (page 12) rendered
and read as an image because the scored bundle's values come from it. Wave B36
(`results/validation/kinetic_core_b36_prereg.md`).

| field | value |
|---|---|
| Title | "The Impact of Varying Enzymatic Pretreatment Durations of Wheat Gluten on the Flavour Characteristics of High-Moisture Plant-Based Extrudates" |
| Venue | Foods 2026, 15, 912 |
| DOI | 10.3390/foods15050912 |
| Group | Zhengzhou University of Light Industry / Henan University of Technology (An, Bai); Massey University (Kaur) |
| Systems | high-moisture extrudates (HMPE) of soy protein isolate : wheat-gluten hydrolysate 6:4 (dry basis); the gluten pre-treated with enzyme for 0, 20, 40, 60 or 80 min. HMPE-0 min is the control and the bundle's system |
| What is measured | free amino acids of the extrudates (Table 1, mg/100 g); HS-SPME-GC-MS volatiles as relative contents against cyclohexanone (Table 2, µg/kg); relative odour activity values (Table 3); sensory and E-nose |
| Bundle | `external_validation_li_2026_spi_wg_hme_control` (external validation, matrix ranking) |

## 1. What this paper is and is not, for this model

A **formulation study at one extrusion condition**: no rate, no barrier, no time course. Its value is
an end-of-process volatile profile of a real high-moisture extrudate, which the bundle scores
externally. Quantification is **semi-quantitative** (one internal standard, cyclohexanone; no
compound-specific calibration), which the bundle's `quantification_class` already records.

## 2. The extrusion, verbatim from sec. 2.3 (page 4)

"co-rotating twin-screw extruder (CLEXTRAL Ev025 ...), featuring a screw diameter of 25 mm and
[L/D 24:1] ... soy protein isolate and wheat gluten hydrolysates at a dry basis ratio of 6:4, was
processed at a screw speed of 280 rpm and a feed rate of 4.6 kg/h, with the mixture's moisture
content maintained at approximately 57%. The extruder's barrel was divided into six zones, with the
temperature profile set at 30, 90, 120, 140, 150, and 160 °C from zone I to VI, respectively. A
cooling die ..." (cooled to 60 °C). The bundle's 160 °C is the last barrel zone. **The paper prints
no residence time and no pH or water activity of the blend**: the bundle's 25 s and its pH 7.0 (the
gluten pretreatment pH, transposed) stay what its own notes call them, assumptions.

## 3. Table 2, HMPE-0 min column (the control) — every row, verbatim, mean ± SD, µg/kg

| compound | RI | threshold µg/kg | HMPE-0 min | HMPE-20 min (for the record) |
|---|---:|---:|---:|---:|
| Hexanal | 780 | 5 | **605.64 ± 6.50** | 522.77 ± 4.30 |
| Heptanal | 902 | 2.8 | 89.88 ± 0.48 | 72.36 ± 0.99 |
| Benzaldehyde | 960 | 750 | 199.14 ± 7.11 | 441.27 ± 4.80 |
| Benzeneacetaldehyde | 1044 | 6.3 | – | 42.39 ± 1.08 |
| Nonanal | 1102 | 1.1 | **74.37 ± 0.11** | 72.66 ± 1.46 |
| Decanal | 1205 | 3 | 29.42 ± 1.08 | 25.28 ± 1.59 |
| 2-Methyl-3-octanone | 985 | 21 | 24.05 ± 1.07 | – |
| 2-Nonanone | 1093 | 41 | 46.47 ± 0.90 | 45.03 ± 2.03 |
| 2-Decanone | 1192 | 8.3 | 34.29 ± 4.77 | 25.80 ± 1.66 |
| 3-Methyl-butanol | 730 | 460 | 14.41 ± 0.20 | 187.82 ± 3.35 |
| 1-Hexanol | 868 | 5.6 | **20.04 ± 0.66** | 35.38 ± 1.58 |
| 1-Octen-3-ol | 978 | 1.5 | 50.23 ± 0.72 | 47.23 ± 0.27 |
| 3,5-Octadien-2-ol | 1037 | – | 29.44 ± 0.63 | 29.38 ± 2.12 |
| 2-Ethylfuran | 691 | 8000 | 168.53 ± 4.81 | 188.98 ± 6.88 |
| Furfural | 836 | 9.56 | **–** (not detected) | 13.15 ± 0.36 |
| 2-Furanmethanol | 852 | 1900 | 231.90 ± 3.41 | 855.95 ± 6.43 |
| 2-Pentylfuran | 994 | 5.8 | **5625.80 ± 63.75** | 5954.76 ± 23.9 |
| Maltol | 1118 | 1.24 | 221.51 ± 5.74 | 353.24 ± 7.12 |
| 3-Hydroxy-2,3-dihydromaltol | 1134 | – | 15.15 ± 0.27 | 62.76 ± 1.77 |
| 3-Phenylfuran | 1228 | – | – | 18.61 ± 1.01 |
| 2-Butylthiophene | 1070 | – | 38.87 ± 0.80 | 38.62 ± 1.47 |
| 2-Pentylthiophene | 1164 | – | 134.60 ± 2.84 | 144.47 ± 2.18 |
| 2-Hexylthiophene | 1274 | – | 30.65 ± 0.75 | 42.66 ± 2.16 |
| Pyrazine | 737 | 2 | 78.57 ± 0.54 | 102.82 ± 4.41 |
| 2,5-Dimethylpyrazine | 916 | 1.75 | 134.30 ± 4.36 | 376.42 ± 10.04 |
| 2-Ethyl-3-methylpyrazine | 1005 | 500 | 60.48 ± 0.49 | 144.04 ± 8.52 |
| 2-Ethenyl-6-methylpyrazine | 1017 | 40 | 11.41 ± 1.72 | 45.02 ± 2.08 |
| 2-Acetylpyrazine | 1022 | 60 | 17.08 ± 0.60 | – |
| 3-Ethyl-2,5-dimethylpyrazine | 1081 | 8.6 | 132.58 ± 3.33 | 174.50 ± 2.29 |
| 2,3-Diethyl-5-methylpyrazine | 1200 | 0.0031 | 18.86 ± 0.12 | 21.87 ± 0.65 |
| n-Pentylpyrazine | 1216 | 1 | 8.17 ± 0.88 | – |
| 2,5-Dimethyl-3-(3-methylbutyl)pyrazine | 1308 | 600 | 34.60 ± 0.72 | 198.75 ± 1.47 |

Bold = the four values the bundle scores.

## 4. The bundle, checked against the print (wave B36)

| bundle value | print, HMPE-0 min | verdict |
|---|---|---|
| hexanal 605.6 | 605.64 ± 6.50 | matches |
| 1-hexanol 20.04 | 20.04 ± 0.66 | matches |
| 2-pentylfuran 5625.8 | 5625.80 ± 63.75 | matches (the 2026-08-27 correction from the maltol row was right) |
| nonanal 72.66 | **74.37 ± 0.11** | **wrong column.** 72.66 ± 1.46 is HMPE-20 min. The 2026-08-27 note fixed the row (decanal → nonanal) and took the column to its right. Corrected by B36 to 74.37. |

**Not added, with the reason.** 2,5-Dimethylpyrazine (134.30 ± 4.36) is a compound the trunk
carries, but the bundle charges only the isolate as a lipid carrier and no free amino acid or sugar,
so the engine's reachability rule would refuse it by name: a refusal is not a validation. Furfural is
not detected in the control and the paper prints no detection limit. Heptanal, benzaldehyde,
1-octen-3-ol, 2-nonanone and the thiophenes are not carried by the trunk. The free amino acids of
the extrudates (Table 1) describe the product after extrusion, not the charge, and are not used.
