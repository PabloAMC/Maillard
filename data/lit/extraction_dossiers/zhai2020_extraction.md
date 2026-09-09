# Zhai et al. 2020 — EXTRACTION (TTCA / ARP proportions of the xylose-cysteine intermediate)

**Source on disk:** `data/articles/zhai2020.pdf` (owner's download, 2026-09-07). Read-only extraction.

| field | value |
|---|---|
| Title | "Transformation between 2-Threityl-thiazolidine-4-carboxylic Acid and Xylose-Cysteine Amadori Rearrangement Product Regulated by pH Adjustment during High-Temperature Instantaneous Dehydration" |
| Authors | Yun Zhai, Heping Cui, Khizar Hayat, Shahzad Hussain, Muhammad Usman Tahir, Shibin Deng, Qiang Zhang, Xiaoming Zhang, Chi-Tang Ho |
| Venue | J. Agric. Food Chem. 2020, 68, 10884-10892 |
| DOI | 10.1021/acs.jafc.0c04287 |

## What it holds

A process paper (spray-drying), not a kinetics paper. The one fact the repo takes: in the xylose-cysteine
"Maillard reaction intermediate" prepared by the group's standard route (Xyl-Cys at pH 7.4, 90 C, then
vacuum dehydration) **TTCA is the predominant form; the Amadori rearrangement product is 6.03 % of the
total without spray-drying**, rising to 20.83 % at a 190 C inlet and to 47.23 % (59.48 % of formation)
when the stock solution is adjusted to pH 9.5 before spray-drying; conversion is favoured over pH
7.5-9.5. The thiazolidine is the stable form under acidic/neutral conditions.

## Consequence for the repo

Wang 2026 prepares its "Cys-Amadori intermediate" by Zhai's route (their Methods cite Zhai 2019) and
characterises it only by mass (C8H15NO6S, m/z 254), which TTCA and the ARP share. Under this paper the
material is ~94 % TTCA, so the Wang pot is charged as **TTCA** (the species the core carries) with the
ARP share recorded as a caveat. Also: the Amadori compound is the more reactive isomer, so a pot with
a higher ARP share browns and degrades faster than a pure TTCA pot at the same loading.
