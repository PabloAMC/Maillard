# Lin, Ting, Ndraha, Hsiao & Sung 2022 — EXTRACTION (chitosan in a fructose–asparagine model solution at 180 °C)

**Source on disk:** `data/articles/lin2022.pdf` (1.7 MB; downloaded 2026-09-11 at the 2026-09-11
reading-list row's request). Read 2026-09-11 via `pdftotext -layout`. Wave B36.

| field | value |
|---|---|
| Title | "Effect of Chitosan Incorporation on the Development of Acrylamide during Maillard Reaction in Fructose–Asparagine Model Solution and the Functional Characteristics of the Resultants" |
| Venue | Polymers 2022, 14, 1565 |
| DOI | 10.3390/polym14081565 |
| Group | National Taiwan Ocean University (Sung) — the 2021 glucose paper's group |
| Systems | 0.5 % fructose + 0.5 % asparagine (± 0.5 % chitosan) in water or 1 % acetic acid, pH 5.8 → 6.0 with NaOH, topped up to 100 mL, "placed in a dry-bath incubator (DB200-2, Yisheng ...)" at 180 °C for 10, 20 and 30 min (sec. 2.2) |
| What is measured | acrylamide (HPLC, calibration 0–3125 ppb), HMF (calibration 0.48–750 ppm), reducing sugar, asparagine, pH, browning, functional properties |
| Bundle | `mp_holdout_fructose_asparagine_180C_Lin2022` (the water arm, 30 min) |

## 1. What this paper is and is not, for this model

The fructose twin of `lin2021_extraction.md`: a **time course in figures**, with the scored
numbers repeated in the prose. The vessel is a tube in a dry-bath incubator; its size and closure
are not printed.

## 2. Numbers printed in the prose (verbatim)

| quantity | value | where |
|---|---|---|
| acrylamide, asparagine–fructose, 30 min, water / 1 % acetic acid | **1859 / 9401 ppb** | Results, "The heating of asparagine-fructose solution for 30 min generated 1859 and 9401 ppb in water and acetic acid, respectively." |
| HMF, asparagine–fructose arm, 30 min | **12.28 ppm** | Results, "(12.28 ppm (Figure 4A) and 5.13 ppm [18])"; 91.56 ppm in the same sentence is the fructose-only arm |

## 3. The bundle, checked against the print (wave B36)

| value | print | verdict |
|---|---|---|
| acrylamide 1859 ppb | 1859 ppb (water) | matches |
| HMF 12 280 ppb | 12.28 ppm | matches |

The charge (27.75 mM fructose, 33.3 mM asparagine) is the bundle's own arithmetic and is unchanged.
