# Ye et al. 2024 — EXTRACTION (tea components against acrylamide in a glucose–asparagine model at 180 °C)

**Source on disk:** `data/articles/ye2024.pdf` (3.6 MB; downloaded 2026-09-11 at the 2026-09-11
reading-list row's request). Read 2026-09-11 via `pdftotext -layout`. Wave B36.

| field | value |
|---|---|
| Title | "Tea's Characteristic Components Eliminate Acrylamide in the Maillard Model System" |
| Venue | Foods 2024, 13, 2836 |
| DOI | 10.3390/foods13172836 |
| Authors | Zhihao Ye, Haojie Xu, Yingying Xie, Ziqi Peng, Hongfang Li, Ruyan Hou, Huimei Cai, Wei Song, Chuanyi Peng, Daxiang Li (Anhui Agricultural University) |
| Systems | "an equimolar solution of glucose and asparagine ... in phosphate buffer (0.1 M, pH 6.86), and 4 mL of the solution was transferred to a 25 mL thick-walled pressurized glass tube in an oil bath ... at 180 °C and kept for 30 min" (sec. 2.2, verbatim in the bundle); then 26 tea components at 0.1–10 g/mol Asn |
| What is measured | acrylamide by LC-QQQ-MS/MS with 13C3-acrylamide, expressed as µmol per mol asparagine (Figure 1 for the control; Table 1 for the 26 components) |
| Bundle | `mp_holdout_glucose_asparagine_180C_Ye2024` (the control, unit `umol_per_mol_limiting_precursor`) |

## 1. What this paper is and is not, for this model

An **inhibitor screen**; the model uses only its control. Two things the print settles: the vessel
sentence the bundle carries is verbatim, and **the reactant molarity is never stated** (the bundle's
0.2 M is taken from the paper's cited method, Knol 2005, and is labelled an assumption in its
`precursor_concentration_provenance`).

## 2. The number printed in the prose (verbatim)

"The level of acrylamide under the control condition (i.e., 0 g/mol Asn) was **140.58 ± 13.92
µmol/mol Asn**" (sec. 3.1, with Figure 1). Inhibition by the tea components ran to 13.22–25.48 % at
the low dose (same paragraph).

## 3. The bundle, checked against the print (wave B36)

140.58 ± 13.92 µmol/mol Asn: matches. Nothing to add — the 26 complexes are inhibitors the model
does not represent.
