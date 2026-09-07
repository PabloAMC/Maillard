# Liu, Wang, Hui, Fang & Zhang 2021 — EXTRACTION (2-furfurylthiol from glucose + cysteine, with ribose)

**Source on disk:** `data/articles/liu2021.pdf` (owner's download, 2026-09-07). Read-only.

| field | value |
|---|---|
| Title | "New insight into the formation mechanism of 2-furfurylthiol in the glucose-cysteine reaction with ribose" |
| Venue | Food Research International 2021, 143, 110295 |
| DOI | 10.1016/j.foodres.2021.110295 |

## 1. Methods

Model 1 (mechanism): glucose 277.00 mg + cysteine 198.00 mg (the Beijing-duck breast levels) in
tripotassium phosphate buffer (5 mL for the duck-level models; the carbon-module models use 77.39 mL
water + a 3.87 mL buffer aliquot), pH 5.88, heated at 168 C for 20 / 40 / 60 min in a heated metal block
under magnetic stirring, cooled in ice. Volatiles by HS-SPME GC-MS, quantified against 2-methyl-3-heptanone
as internal standard (semi-quantitative, "ng" per vial). Six duck-level models (Cys; Glc; Rib; Glc+Cys;
Rib+Cys; Glc+Rib+Cys) at 168 C / 60 min.

## 2. Table 3 — glucose-cysteine, ng per vial (mean +/- SD, n = 3)

| compound | 20 min | 40 min | 60 min |
|---|---:|---:|---:|
| 2-furfurylthiol | 44.73 +/- 10.00 | 107.15 +/- 9.59 | 133.43 +/- 17.30 |
| 2-furfural | 0 | 4.23 +/- 0.40 | 4.09 +/- 0.77 |
| 2-furanmethanol | 3.73 +/- 0.38 | 35.40 +/- 5.45 | 81.65 +/- 2.07 |

Carbon-module labelling: FFT's furan ring comes from glucose fragmentation (unlabeled [M]+ 42-54 % at
20 min), with 2-furfural and 2-furanmethanol as intermediates and 2-furanmethanol + H2S the last step.
Adding ribose at the duck level raises FFT sharply (Table 6; the ribose route dominates).

## 3. What the repo takes

A HEXOSE-only FFT time course (glucose + cysteine, no pentose) — the entry the sulfur lane declares
UNIDENTIFIED (hexose-only charges predict zero MFT since B9 and refuse the thiol targets). Directional
claim HEX-T-01 (FFT rises 20 -> 60 min at 168 C) is recorded NOT EVALUABLE by the engine's own declaration;
it becomes the first validation row for whichever wave installs a hexose -> furfural/furfuryl alcohol
-> FFT route. Levels are IS-relative and not transferable.
