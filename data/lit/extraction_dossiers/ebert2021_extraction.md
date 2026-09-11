# Ebert, Michel, Nedele, Baune, Terjung, Zhang, Gibis & Weiss 2021 — EXTRACTION (two pea protein isolates and their dry and wet texturates, as normalized peak-area percentages)

**Source on disk:** `data/articles/ebert2021.pdf` (1.2 MB; downloaded 2026-09-11). Read 2026-09-11
via `pdftotext -layout` with Tables 1 and 3 verified against page renders. Wave B37.

| field | value |
|---|---|
| Title | "Influence of protein extraction and texturization on odor-active compounds of pea proteins" |
| Venue | Journal of the Science of Food and Agriculture; the PDF on disk is the accepted version and prints **no volume and no page numbers**, only "J Sci Food Agric 2021" |
| DOI | 10.1002/jsfa.11437 |
| Group | University of Hohenheim with the German Institute for Food Technology (DIL) |
| Samples | Pea Protein I (Pisane C9) and II (M9), each with a **dry texturate (TVP)** and a **wet texturate (WTP)** made by DIL |
| Quantification | direct-immersion stir-bar sorptive extraction (PDMS Twister), 5 mL of a 1 % dispersion, 1000 rpm, 2 h at 25 °C; **normalized peak-area percentages (the 100 % method) over a 0–38 min window** |

## 1. Why it was fetched, and why it cannot do the job

It was fetched because it prints a pea protein isolate beside its own extrudates — the unheated
column the experiments guide asks for. It does. **But two things disqualify it as a benchmark.**

1. **The values are not concentrations.** There is no internal standard, no calibration curve and no
   isotope dilution; every number is that compound's share of the total integrated chromatogram
   area. A value can move because the compound moved or because the denominator moved. The authors
   say so and call for absolute work.
2. **The paper prints no process parameters at all.** No extrusion temperature, no residence time, no
   moisture of the premix, no screw speed — the texturates were made by a third party and are
   described only as "dry texturization" and "wet texturization". There is no time–temperature
   coordinate anywhere in the paper, so no pot can be charged from it.

## 2. Table 3 — peak area (%), the rows the trunk carries

| compound | Pea Protein I | Pea TVP I | Pea WTP I | Pea Protein II | Pea TVP II | Pea WTP II |
|---|---:|---:|---:|---:|---:|---:|
| **Hexanal** | 3.29 ± 1.05 | 3.16 ± 0.03 | 0.52 ± 0.02 | 4.40 ± 0.33 | 2.37 ± 0.06 | 3.50 ± 0.13 |
| Octanal | 0.88 ± 0.27 | 1.04 ± 0.05 | n.i. | 0.76 ± 0.01 | 0.78 ± 0.02 | 1.32 ± 0.19 |
| **Nonanal** | 2.36 ± 0.58 | 2.20 ± 0.11 | 0.58 ± 0.07 | 1.69 ± 0.14 | 1.88 ± 0.18 | 2.57 ± 0.05 |
| Decanal | 0.59 ± 0.26 | 0.98 ± 0.02 | n.i. | n.d. | n.d. | n.d. |
| (E,E)-2,4-Decadienal | 4.71 ± 3.54 | 2.61 ± 0.37 | 3.15 ± 1.07 | 5.15 ± 0.38 | 0.94 ± 0.11 | 7.85 ± 0.04 |
| **2-Pentylfuran** | 5.42 ± 0.61 | 3.28 ± 0.46 | 5.95 ± 1.08 | 5.96 ± 0.08 | 1.45 ± 0.03 | 7.30 ± 0.47 |
| (E,E)-3,5-Octadien-2-one | 4.36 ± 0.12 | 2.89 ± 0.11 | 4.80 ± 0.01 | 5.72 ± 0.35 | 1.98 ± 0.07 | 4.11 ± 0.01 |
| total peak area (×10⁶) | 247.49 ± 7.93 | 187.16 ± 1.40 | 216.72 ± 33.45 | 241.14 ± 40.38 | 183.82 ± 11.43 | 249.95 ± 9.29 |

`n.i.` = not integrable (below 65 000 per minute); `n.d.` = not detectable in the total ion current.
**These are different facts and are not collapsed to zero.** Note that **Table 1 and Table 3 use
different column orders** — the single likeliest transcription error in this paper.

Composition (Table 1, manufacturer's figures, no SDs): Pea Protein I / II are 6.0 / 5.7 % moisture,
76.4 / 77.1 % crude protein, 9.0 / 8.0 % total fat of which 5.4 / 4.6 % polyunsaturated. The wet
texturates are 63.8 / 64.6 % moisture.

**Absent from the paper entirely**: 1-hexanol, heptanal, benzaldehyde, 1-octen-3-ol,
2-acetyl-1-pyrroline, methional, furfural, furaneol. The pyrazines — 2,5-dimethylpyrazine,
3-ethyl-2,5-dimethylpyrazine, 2-butyl-3,5-dimethylpyrazine — are **identified but never quantified**;
they are stars on Figure 2. That is the one Maillard signal in the paper and it has no number.

## 3. What the repository takes

A **direction**, recorded and not scored: wet texturization cuts hexanal roughly sixfold in Pea
Protein I (3.29 → 0.52 % of area) and 2-pentylfuran by 40 % (TVP I) and 75 % (TVP II), while dry
texturization is where pyrazines appear. Nothing else. The paper's own kinetic content is a citation
to Zamora, Navarro, Aguilar & Hidalgo, "Lipid-derived aldehyde degradation under thermal conditions",
Food Chem. 174:89–96 (2015) — **that** is the lead for a thermal aldehyde-degradation law, and it is
not on disk.
