# Bischof et al. 2023 — EXTRACTION (beef loin, 1H NMR, glucose and nucleotides over 28 days)

**Source on disk:** none. Read 2026-09-14 from the publisher's HTML full text
(`link.springer.com/article/10.1007/s00217-023-04283-0`) through an automated fetch. The read
returned Table 1 values but repeated the same numbers under the dry-aged and wet-aged headings,
which is almost certainly the reader collapsing columns; only the day-0 and day-28 means below are
carried and they must be checked against the PDF before any panel use. Written for the
cultivated-tissue composition box, not for the core fit.

| field | value |
|---|---|
| Title | "NMR-based comparison of the metabolome of beef from Simmental and black-and-white young bulls during wet- and dry-aging" |
| Venue | European Food Research and Technology 2023, 249, 2113-2124 |
| DOI | 10.1007/s00217-023-04283-0 |

## 1. Methods

M. longissimus thoracis et lumborum from 7 Simmental and 7 Black-and-White young bulls, dry- and
wet-aged, sampled at 0, 7, 14, 21 and 28 days. Polar extract, 1H NMR, 30 metabolites quantified.
Table 1: "Concentration of metabolites (mean ± SD in µmol/g wet meat)".

## 2. Findings that matter (µmol per g wet meat, mean ± SD)

| compound | Simmental d0 | Simmental d28 | Black-and-White d0 | Black-and-White d28 |
|---|---|---|---|---|
| glucose | 1.83 ± 0.91 | 2.36 ± 1.15 | 4.16 ± 1.48 | 4.05 ± 1.59 |
| IMP | 5.07 ± 1.05 | 1.84 ± 0.57 | 5.36 ± 0.99 | 2.31 ± 0.88 |
| inosine | 0.29 ± 0.12 | 0.73 ± 0.27 | 0.39 ± 0.13 | 0.93 ± 0.36 |
| hypoxanthine | 0.39 ± 0.15 | 1.76 ± 0.53 | 0.41 ± 0.14 | 1.58 ± 0.57 |
| leucine | 0.68 ± 0.13 | 1.47 ± 0.31 | 0.59 ± 0.11 | 1.39 ± 0.26 |
| methionine | 0.36 ± 0.08 | 0.53 ± 0.11 | 0.31 ± 0.07 | 0.48 ± 0.09 |

Ribose is not among the 30 metabolites quantified.

## 3. What the repo takes

Beef-side ranges for the composition box, in mM in tissue water (µmol/g ÷ 0.75), spanning
mean − SD to mean + SD across both breeds and both ages:

- **glucose (free)**: 1.2 to 7.5 mM.
- **leucine (free)**: 0.64 to 2.4 mM (agrees with Muroya 2019's 0.35 to 1.1 mM at the low end).
- **IMP**: 1.7 to 8.5 mM (inside Muroya 2019's 0.10 to 10 mM).

Fourteen animals, two breeds, one laboratory. The day-0 glucose difference between breeds (1.8 vs
4.2 µmol/g) is larger than the ageing change within a breed.
