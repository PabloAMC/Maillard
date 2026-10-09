# Dashmaa et al. 2026 — EXTRACTION (review: flavoromics of meats; carries Aliani 2013's sugar values, which are CHICKEN)

Stem `hwang2026` after the corresponding (last) author, Inho Hwang; the first author is Dashdorj Dashmaa and
the PDF on disk is named after her. `scripts/reading_audit.py` carries the alias.

**Source on disk:** `data/articles/Dashmaa2026.pdf` (the publisher's open-access PDF). First read 2026-09-14
from the publisher's HTML through an automated fetch; the one sentence the repo took, and the reference it
cites, checked against the PDF by eye the same day (p. 6 of the PDF, the paragraph on sugar-related
metabolites; reference list). A REVIEW, used here only as the carrier of one primary number set that the repo has not read at
source. Written for the cultivated-tissue composition box, not for the core fit.

| field | value |
|---|---|
| Title | "Flavoromics of meats as a function of postmortem proteolysis" |
| Authors | D. Dashmaa, L. Xi, H. Van Ba, J. Dawoon, I. Hwang |
| Venue | Food Science of Animal Resources 2026, 46:46 |
| DOI | 10.1007/s44463-025-00015-y |

## 1. What it carries

"Naturally, beef contains approximately 1.48, 0.33, and 0.26 mg/g of glucose, fructose, and
ribose, respectively", attributed to Aliani et al. 2013. **The PDF's reference list gives that paper as: Aliani, M., Farmer, L. J.,
Kennedy, J. T., Moss, B. W., & Gordon, A. (2013). Post-slaughter changes in ATP metabolites, reducing and
phosphorylated sugars in CHICKEN meat. Meat Science, 94, 55-62.** The review's sentence says "beef"; its
source measured chicken. The muscle, ageing state and moisture basis are unknown here, and the species
attribution is the review's error.

Also: "The ATP concentration is usually approximately 5-8 µmol/g in resting muscle" (Feiner 2006),
consistent with Muroya 2019's D0 ATP of 6.5 µmol/g.

## 2. What the repo took on 2026-09-14 as SECONDARY values, and no longer takes (same day, PDF check)

Converted to mM in tissue water at 75 % moisture (mg/g ÷ molar mass × 1000 ÷ 0.75):

- **ribose, beef**: 0.26 mg/g ÷ 150.1 → 1.73 mmol/kg → 2.3 mM. The upper corner of the box's beef
  ribose range; the lower corner is one sixth of it, from Koutsidis 2008b's abstract
  (`koutsidis2008b_extraction.md`).
- **glucose, beef**: 1.48 mg/g ÷ 180.2 → 8.2 mmol/kg → 11 mM. Above Bischof 2023's primary
  1.2 to 7.5 mM; carried as the upper corner of the glucose range so the box spans both.

**Superseded the same day.** Once Koutsidis 2008a and 2008b were read from their PDFs
(`koutsidis2008a_extraction.md`, `koutsidis2008b_extraction.md`) and Bischof 2023's glucose rows were
corrected (`bischof2023_extraction.md`), every beef sugar range in the box rests on a primary table and
this review is cited on none of them. The two cited points survive only as cross-checks: 2.3 mM ribose
sits at the top edge of the primary 0.33 to 2.2 mM, and 11 mM glucose inside the primary 2.4 to 15 mM.
Aliani 2013 itself has no PDF on disk and no dossier, and needs none for the beef box: the PDF's reference
list shows it is a chicken paper, so the two cited points were never beef values and are not even
cross-checks. That the chicken numbers happen to fall inside the beef ranges is a coincidence of similar
muscle chemistry, not a confirmation.
