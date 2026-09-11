# Kong, Wu, Li, Kang, Wang, Xie & Yu 2024 — EXTRACTION (soy protein isolate, one 120 °C thermal point against an untreated control, with ultrasound variants)

**Source on disk:** `data/articles/Kong2024.pdf` (3.1 MB; downloaded 2026-09-11). Read 2026-09-11 via
`pdftotext -layout`, Table 1 checked against a rendered page image. Wave B37.

| field | value |
|---|---|
| Title | "Analyzing changes in volatile flavor compounds of soy protein isolate during ultrasonic-thermal synergistic treatments using electronic nose and HS-SPME-GC-MS combined with chemometrics" |
| Venue | Food Chemistry 445 (2024), article 138795 |
| DOI | 10.1016/j.foodchem.2024.138795 |
| Group | Northeast Agricultural University, Harbin |
| Samples | **Control** (untreated); **TT** = oil bath 120 °C for 150 s; **UT-250 / UT-350 / UT-450** = 20 kHz ultrasound at 250/350/450 W for 7 min, then 120 °C for 150 s |
| Quantification | HS-SPME (50/30 µm DVB/CAR/PDMS), 1.0 g powder, 60 °C, 40 min; **2-methyl-3-heptanone internal standard with NIST 14L library identification and no calibration curves** |

## 1. What class of number this is

The units are printed as µg/kg, but there is **no authentic-standard calibration and no response
factor** — each value is the analyte peak area scaled to one unrelated internal standard.
These are **internal-standard equivalents**, one tier below `xu2024_extraction.md` and one tier above
a bare peak area. The repository's `quantification_class` vocabulary already has the right label for
this: semi-quantitative against an internal standard.

## 2. Table 1 — contents (µg/kg), the requested rows, verbatim

| # | compound | Control | TT (120 °C, 150 s) | UT-250 | UT-350 | UT-450 |
|---:|---|---:|---:|---:|---:|---:|
| 1 | **Hexanal** | 177.72 ± 2.01 | 134.84 ± 3.02 | 85.53 ± 1.11 | 52.35 ± 2.22 | 74.40 ± 1.42 |
| 2 | (E)-2-Hexenal | 75.26 ± 5.23 | 13.19 ± 1.04 | 4.89 ± 0.94 | 3.31 ± 0.34 | 3.32 ± 0.42 |
| 3 | **Heptanal** | 8.00 ± 1.03 | 13.74 ± 2.33 | 9.0 ± 1.05 | 5.91 ± 0.87 | 6.09 ± 0.95 |
| 4 | (Z)-2-Heptenal | 4.77 ± 1.03 | 3.21 ± 1.03 | 1.52 ± 0.23 | – | – |
| 5 | **Benzaldehyde** | 12.15 ± 1.09 | 12.08 ± 1.18 | 7.84 ± 0.12 | 5.87 ± 0.23 | 6.45 ± 0.11 |
| 6 | **Octanal** | 11.66 ± 1.22 | 6.66 ± 1.18 | 20.80 ± 3.32 | 21.02 ± 2.27 | 18.14 ± 2.11 |
| 7 | **Nonanal** | 164.58 ± 4.66 | 94.97 ± 3.23 | 145.04 ± 5.12 | 135.90 ± 5.77 | 126.02 ± 4.55 |
| 8 | Pentanal | – | – | 8.80 ± 0.71 | – | – |
| 9 | (E)-2-Nonenal | 5.49 ± 0.41 | – | 4.47 ± 0.33 | – | – |
| 10 | **Decanal** | 62.89 ± 1.27 | 34.86 ± 0.72 | 63.75 ± 4.29 | 75.79 ± 4.44 | 73.37 ± 3.41 |
| 12 | (E)-2-Octenal | 3.67 ± 0.01 | 4.19 ± 0.01 | 5.31 ± 0.02 | 5.69 ± 0.04 | 5.33 ± 0.01 |
| 16 | 1-Pentanol | 6.23 ± 1.11 | 5.88 ± 2.01 | 4.17 ± 0.31 | 3.30 ± 0.43 | 3.33 ± 0.54 |
| 17 | **1-Octen-3-ol** | 14.47 ± 1.44 | 11.06 ± 2.31 | 7.30 ± 0.31 | 5.61 ± 0.44 | 6.00 ± 0.32 |
| 33 | **2-Pentylfuran** | 10.28 ± 0.19 | 9.88 ± 0.22 | – | – | – |

"–" is the paper's "not recognized", which is **not** a measured zero. **1-Hexanol is not in this
paper**, and there is no pyrazine anywhere in it. The heptenal here is the **(Z)-2** isomer, not
Xu 2024's (E)-2 — the two must never be merged.

## 3. What it supports and what it does not

Same direction as Xu 2024 and independently: **heating a soy protein isolate lowers its lipid-derived
volatiles** — hexanal 177.72 → 134.84, nonanal 164.58 → 94.97, decanal 62.89 → 34.86 µg/kg at
120 °C for 150 s. The abstract puts it verbatim for the ultrasound-plus-heat arm: hexanal,
(E)-2-hexenal and 1-octen-3-ol "reduced by 70.60 %, 95.60 % and 61.23 %".

**One thermal point only** (120 °C, 150 s): no time course, no temperature series, no kinetics. The
three UT columns confound heat with seven minutes of 20 kHz cavitation and are not thermal points.

**An unresolved confound that blocks benchmark use.** The methods say every sample was vacuum
degassed for 10 minutes at 0.1 MPa and spray dried, then name only the four treated samples; whether
the control received the degassing is not stated, and the supplementary flowchart that would settle
it is not in the PDF on disk. A ten-minute vacuum strips volatiles on its own, so a control-versus-TT
difference cannot be attributed to the heat until that is resolved. Recorded, not resolved.

Table 2 (odour activity values and contribution rates) is transcribed in the repository's reading
notes but is derived from Table 1 and a threshold set that differs from Xu 2024's; the two threshold
sets must not be mixed.
