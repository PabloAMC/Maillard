# Lombardi-Boccia, Lanzi & Aguzzi 2005 — EXTRACTION (thiamine, riboflavin, niacin and trace elements in raw and cooked meat cuts)

**Source on disk:** `data/articles/lombardi-boccia2005.pdf` (the publisher's PDF). Table 2 read from the PDF
and checked by eye on 2026-09-14. This replaces the excerpt-only version written earlier the same day from
two sentences on a figure-sharing page; the excerpt's two numbers (0.01 and 0.08 mg/100 g) were correct,
the title, venue and DOI it left unconfirmed are confirmed below, and the per-cut table is now on file.
Written for the cultivated-tissue composition box, not for the core fit.

| field | value |
|---|---|
| Title | "Aspects of meat quality: trace elements and B vitamins in raw and cooked meats" |
| Authors | G. Lombardi-Boccia, S. Lanzi, A. Aguzzi (INRAN, Rome) |
| Venue | Journal of Food Composition and Analysis 2005, 18(1), 39-46 |
| DOI | 10.1016/j.jfca.2003.10.007 |

## 1. Methods

Retail and producer meat cuts from Italy: beef (sirloin, fillet, roast beef, topside, thick flank), veal,
lamb, horse, ostrich, pork, chicken, turkey, rabbit. Raw and cooked (pan or grill to the disappearance of
red colour). Thiamine, riboflavin and niacin by HPLC after acid hydrolysis (Barna & Dworschák 1994), so
the thiamine figure is TOTAL thiamine (free plus phosphorylated), which is the pool the sulfur lane's
thiamine route starts from. Table 2 prints mean ± SD in mg per 100 g fresh weight, raw and cooked, with a
weight-loss column; letters mark p < 0.05 among cuts of a species.

## 2. Findings that matter (Table 2, thiamine, mg/100 g fresh weight)

| beef cut | raw | cooked | weight loss on cooking |
|---|---|---|---|
| sirloin | 0.02 ± 0.01 | nd | 37.7 % |
| fillet | 0.08 ± 0.01 | nd | 38.2 % |
| roast beef | 0.05 ± 0.01 | nd | 39.2 % |
| topside | 0.08 ± 0.01 | nd | 43.4 % |
| thick flank | 0.01 ± 0.01 | nd | 40.0 % |

For comparison: veal fillet 0.11, pork cuts 0.6 to 0.9, horse 0.18, lamb 0.16 mg/100 g raw. Riboflavin in
beef 0.09 to 0.17, niacin 5.0 to 6.5 mg/100 g.

Two things the excerpt did not show. First, beef thiamine after cooking was **not detected in any of the
five cuts** (veal likewise), while riboflavin and niacin survived at 40 to 80 %: at these low starting
levels the cooking step consumed the whole pool, which is the behaviour the thiamine route in the sulfur
lane presupposes and a qualitative check on it. Second, the spread across cuts (eightfold) is larger than
any ageing or breed effect the other beef sources print for any precursor.

## 3. What the repo takes

Beef total thiamine 0.01 to 0.08 mg/100 g across five raw cuts, one laboratory. In mM in tissue water at
75 % moisture (mg/100 g × 10 ÷ 300.8 g/mol ÷ 0.75), spanning the lowest mean to the highest mean + SD:
**0.00044 to 0.0040 mM**, labelled `sourced`. The earlier `secondary` range (0.00044 to 0.0049 mM) had its
upper corner from a review (Lee, Lee & Jo 2025, Food Sci Anim Resour 45:303, citing Ramalingam 2019 at
0.08 to 0.11 mg/100 g); that review is unread and no longer cited on the range. Its top sits 20 % above the
primary corner and is not a disagreement.

Nothing here for thiamine in cultured cells. DMEM carries about 4 mg/L thiamine hydrochloride (≈ 12 µM);
what a washed, differentiated construct retains is unmeasured in the public record as of this read.
