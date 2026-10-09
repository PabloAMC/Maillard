# Muroya et al. 2019 — EXTRACTION (beef loin, CE-TOFMS, absolute pools at 0, 1 and 14 days)

**Source on disk:** `data/articles/Muroya2019.pdf` (the publisher's open-access PDF). First read
2026-09-14 through an automated fetch of the same PDF; Table 1 then checked against the PDF by eye
later the same day: **every number in section 2 matched** (means and SEs, all eleven rows). Two
additions from the check are marked below. Written for the cultivated-tissue composition box
(`src/cultivated_tissue_invariance.py`), not for the core fit.

| field | value |
|---|---|
| Title | "Metabolomic approach to key metabolites characterizing postmortem aged loin muscle of Japanese Black (Wagyu) cattle" |
| Venue | Asian-Australasian Journal of Animal Sciences 2019, 32(8), 1172-1185 |
| DOI | 10.5713/ajas.18.0648 |

## 1. Methods

Longissimus thoracis from n = 3 Japanese Black steers (28 months, 632-739 kg), stored at 2 °C, sampled
at 0 (30 min), 1 and 14 days post-mortem. **Added at the PDF check:** the muscle blocks were 36.5 %
crude fat, 14.2 % protein and 47.7 % moisture (Wagyu marbling); intramuscular fat was removed by
hand before the lean pieces were frozen, and the nmol/g are per g of that trimmed lean, whose own
moisture is not printed. The 75 % moisture used below is the usual lean-beef figure; if the trimmed
lean still carried fat, the mM in tissue water would be higher, by a factor well under two. Capillary-electrophoresis time-of-flight MS with absolute quantification against
standards. Table 1 prints "Mean (nmol/g)" with SE; compounds not detected at any time point are
omitted from the table.

## 2. Findings that matter (Table 1, nmol per g wet tissue, mean ± SE)

| compound | D0 | D1 | D14 |
|---|---|---|---|
| glucose 6-phosphate | 10 707 ± 826 | 16 657 ± 3 995 | 8 777 ± 1 545 |
| ribose 5-phosphate | 0 | 56.6 ± 8.6 | 70.0 ± 16.6 |
| ribulose 5-phosphate | 34.9 ± 8.9 | 180 ± 6 | 215 ± 49 |
| IMP | 78.4 ± 39.2 | 7 574 ± 1 402 | 3 438 ± 206 |
| inosine | 4.5 ± 2.5 | 544 ± 92 | 1 109 ± 240 |
| hypoxanthine | 9.9 ± 2.3 | 335 ± 37 | 2 213 ± 203 |
| ATP | 6 535 ± 452 | 21.1 ± 7.0 | 16.4 ± 2.7 |
| cysteine | 1.6 ± 0.8 | 23.6 ± 6.9 | 107 ± 19 |
| methionine | 44.0 ± 4.1 | 42.6 ± 2.4 | 324 ± 71 |
| leucine | 263 ± 26 | 315 ± 19 | 827 ± 133 |
| glycine | 1 083 ± 16 | 1 121 ± 165 | 1 273 ± 63 |

Free glucose, free ribose and thiamine are not in Table 1 (confirmed at the PDF check). Thiamine
appears only in Table 3 and Figure 4 as a RELATIVE content (peak area, no standard) that rose over
aging; Table 3 also prints relative contents for Cys, Leu and Met that are not the absolute values
and must not be mixed with Table 1.

## 3. What the repo takes

Beef-side ranges for the composition box, converted to mM in tissue water at 75 % moisture
(nmol/g ÷ 1000 ÷ 0.75):

- **cysteine (free)**: 0.002 (D0) to 0.14 (D14) mM. The D0 value is near the detection floor and
  the SE is half the mean; the range is the ageing span, not a laboratory spread.
- **leucine (free)**: 0.35 to 1.1 mM.
- **IMP**: 0.10 (D0, pre-rigor) to 10 mM (D1 peak); 4.6 mM at D14.
- **ribose 5-phosphate**: 0 (D0) to 0.093 mM; the box's lower corner is set at 0.005 mM because the
  draw is log-uniform and D0 is a non-detect, not a zero.

Three animals, one breed, one laboratory. Every range above is a span over ageing time in three
steers, and is used as such: it brackets what a beef reference could be, it does not estimate a
population.
