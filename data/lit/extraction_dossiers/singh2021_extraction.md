# Singh, Shi, Magreault, Kitts, Jarzębski, Siejak & Pratap-Singh 2021 — EXTRACTION (rapid HS-SPME-GC-MS odour-activity screen of soy, pea and brown-rice protein powders)

**Source on disk:** `data/articles/singh2021.pdf` (1.4 MB; downloaded 2026-09-11 at the
2026-09-11 reading-list row's request). Read 2026-09-11 via `pdftotext`; Table 1 (page 3)
rendered and read as an image because its text layer is doubled. Wave B36.

| field | value |
|---|---|
| Title | "A Rapid Gas-Chromatography/Mass-Spectrometry Technique for Determining Odour Activity Values of Volatile Compounds in Plant Proteins: Soy, and Allergen-Free Pea and Brown Rice Protein" |
| Venue | Molecules 2021, 26, 4104 |
| DOI | 10.3390/molecules26134104 |
| Group | University of British Columbia (Kitts, Pratap-Singh); Poznań University of Life Sciences |
| Systems | three commercial protein powders, **never heated**: 1 g powder in 7 mL water, HS-SPME (DVB/CAR/PDMS) "under agitation for 10 min at 40 °C", GC-MS; hexanal-d12 internal standard |
| What is measured | Table 1: three fully quantified standards (hexanal, 2-nonanone, hexanol) in ppb; everything else semi-quantified as equivalents of those standards (Table S1, with odour thresholds and OAVs) |
| Bundles | `pea_isolate_40C_PratapSingh2021`, `soy_isolate_40C_PratapSingh2021` (matrix headspace; **no cook**, refused since B31 unless a starting state is declared) |

## 1. What this paper is and is not, for this model

A **method paper with a three-powder survey**. There is no thermal process: the bundles' 40 °C /
10 min is the SPME incubation, as their vessel notes have said since 2026-09-04, and the B31 rule
refuses to answer a formation from zero for a pot that was never cooked. The paper's numbers are
therefore **starting levels of a raw isolate**, the very thing clause 2 of the B31 refusal asks a
user to declare. Its remaining use is as that declaration for a pea or soy isolate of this kind.

## 2. Table 1 — fully quantified (ppb, mean ± SD), verbatim from the image

| compound | Pea | Brown rice | Soy |
|---|---:|---:|---:|
| Hexanal | **1138.00 ± 297.30** | 22,590.24 ± 1643.70 | **1621.71 ± 159.69** |
| 2-Nonanone | 6.382 ± 0.62 | 94.02 ± 12.38 | n.d. |
| Hexanol | n.d. | 102.04 ± 9.30 | n.d. |

(The text layer prints "1643.70" next to the soy hexanal; the image shows it is the brown-rice SD.)

## 3. Semi-quantified numbers printed in the prose (hexanal equivalents, ppb)

| quantity | value | where |
|---|---|---|
| 2-pentylfuran, soy | **2492 ± 199** | sec. 2.4 ("2-pentyl furan was found in large quantities in soy and pea proteins (2492 ± 199 and 638 ± 49 ppb equivalents of hexanal, respectively ...)") |
| 2-pentylfuran, pea | **638 ± 49** | same sentence |
| 2-n-butylfuran, brown rice | 386 ± 30 | same paragraph |
| total aldehydes: brown rice / pea / soy | 40,250 ± 3938 / 1359 ± 321 / 1998 ± 201 | sec. 2.2 |
| soy alcohol fraction | 40 ± 9 of 1-octen-3-ol only | sec. 2.3 (as the soy bundle's correction note already quotes) |

## 4. The bundles, checked against the print (wave B36)

| bundle value | print | verdict |
|---|---|---|
| pea hexanal 1138.0 ± 26 % | 1138.00 ± 297.30 (26.1 %) | matches |
| soy hexanal 1621.71 ± 10 % | 1621.71 ± 159.69 (9.8 %) | matches |
| pea 2-pentylfuran 638 ± 8 % | 638 ± 49 (7.7 %), **hexanal-equivalent, semi-quantified** | matches; the semi-quantitative basis is now stated in the vessel note |
| soy 2-pentylfuran 2492 ± 8 % | 2492 ± 199 (8.0 %), hexanal-equivalent | matches |

Nothing to add: hexanol is n.d. in both scored powders (the 2026-08-27 removals were right);
2-nonanone and 1-octen-3-ol are not carried by the trunk.
