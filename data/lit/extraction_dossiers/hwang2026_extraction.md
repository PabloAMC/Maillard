# Hwang et al. 2026 — EXTRACTION (review: flavoromics of meats; carries Aliani 2013's beef sugar values)

**Source on disk:** none. Read 2026-09-14 from the publisher's HTML through an automated fetch.
A REVIEW, used here only as the carrier of one primary number set that the repo has not read at
source. Written for the cultivated-tissue composition box, not for the core fit.

| field | value |
|---|---|
| Title | "Flavoromics of meats as a function of postmortem proteolysis" |
| Venue | Food Science of Animal Resources, 2026 |
| DOI | 10.1007/s44463-025-00015-y |

## 1. What it carries

"Naturally, beef contains approximately 1.48, 0.33, and 0.26 mg/g of glucose, fructose, and
ribose, respectively", attributed to Aliani et al. 2013. The primary paper has not been read; the
muscle, ageing state and moisture basis of those numbers are unknown here.

Also: "The ATP concentration is usually approximately 5-8 µmol/g in resting muscle" (Feiner 2006),
consistent with Muroya 2019's D0 ATP of 6.5 µmol/g.

## 2. What the repo takes, as SECONDARY values

Converted to mM in tissue water at 75 % moisture (mg/g ÷ molar mass × 1000 ÷ 0.75):

- **ribose, beef**: 0.26 mg/g ÷ 150.1 → 1.73 mmol/kg → 2.3 mM. The upper corner of the box's beef
  ribose range; the lower corner is one sixth of it, from Koutsidis 2008b's abstract
  (`koutsidis2008b_extraction.md`).
- **glucose, beef**: 1.48 mg/g ÷ 180.2 → 8.2 mmol/kg → 11 mM. Above Bischof 2023's primary
  1.2 to 7.5 mM; carried as the upper corner of the glucose range so the box spans both.

Both are labelled `secondary` in the box: a number read from a review, not from the paper that
measured it.
