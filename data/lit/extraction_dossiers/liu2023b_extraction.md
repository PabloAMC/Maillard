# Liu et al. 2023 (LWT) — EXTRACTION (cysteine + ribose at 168 C, 20-60 min, carbon-module labelling)

**Source on disk:** `data/articles/Liu2023b.pdf` (owner's download, 2026-09-07; `Liu2023.pdf` is the
different Food Chem 406:134998 paper). Read-only.

| field | value |
|---|---|
| Title | "Carbon module labeling approach combining with gas chromatography-olfactometry-mass spectrometry technology clarifies the formation mechanism of furan and thiophene derivatives: Ribose and cysteine as a classic case study" |
| Venue | LWT 2023, 182, 114874 |
| DOI | 10.1016/j.lwt.2023.114874 |

## 1. Methods

Ribose + cysteine at the Beijing-duck weight ratio (ribose : cysteine 1 : 10.10; 19.6 mg ribose to
198 mg cysteine per model), in 0.5 mol/L tripotassium phosphate buffer, pH 5.88, glass vial, 168 C for
20 / 40 / 60 min, cooled in cold water. Volatiles by HS-SPME GC-O-MS with 2-methyl-3-heptanone as
internal standard ("ng" per vial, semi-quantitative). [13C5]-ribose 1:1 for the labelling models.

## 2. Table 2 — ng per vial (mean +/- SD)

| compound | 20 min | 40 min | 60 min |
|---|---:|---:|---:|
| 2-methyl-3-furanthiol | 1252.67 +/- 40.70 | 760.33 +/- 17.21 | 628.00 +/- 70.34 |
| 2-furfurylthiol | 175.36 +/- 15.40 | 197.68 +/- 6.64 | 169.51 +/- 10.62 |
| 2-furfural | 5.78 +/- 0.13 | 3.66 +/- 0.96 | 4.06 +/- 0.75 |

Labelling: MFT is mostly the intact C5 ribose skeleton (~81 %) with ~19 % from recombined fragments;
4-hydroxy-5-methyl-3(2H)-furanone is its intermediate; 3-thiophenethiol arises by aldol condensation.

## 3. What the repo takes

A SINK-shape claim at high temperature: at 168 C MFT FALLS by half between 20 and 60 min while FFT is
flat (peak at 40 min within SD) — the first measurement in the corpus of MFT declining with time in a
cysteine-ribose pot, and a direct test of the thiol-sink barrier the fit holds at its 102 kJ/mol ceiling.
Directional claims RIB-T-01 (MFT decreasing 20 -> 60 min) and RIB-T-02 (FFT flat over the same window),
independent (Beijing lab, not in any fit). The pot is cysteine-rich (ribose : cysteine ~ 1 : 15 molar),
so the sugar is the limiting charge.
