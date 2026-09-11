# Lin, Chan, Kao & Sung 2021 — EXTRACTION (HMF and low-molecular-weight chitosan in glucose–asparagine model systems at 180 °C; the paper the three "Chang 2021" bundles score)

**Source on disk:** `data/articles/lin2021.pdf` (2.0 MB; downloaded 2026-09-11 at the 2026-09-11
reading-list row's request). Read 2026-09-11 via `pdftotext -layout`. Wave B36.

| field | value |
|---|---|
| Title | "Effect of Hydroxymethylfurfural and Low-Molecular-Weight Chitosan on Formation of Acrylamide and Hydroxymethylfurfural during Maillard Reaction in Glucose and Asparagine Model Systems" |
| Venue | Polymers 2021, 13, 1901 |
| DOI | 10.3390/polym13121901 |
| Authors | Hong-Ting Victor Lin, Der-Sheng Chan, Ling-Yu Kao, Wen-Chieh Sung (National Taiwan Ocean University). **"Chang" in the bundle names is Chang et al., the method reference [16] the paper follows; the bundles' `source_doi` is this paper and was always right.** |
| Systems | 0.5 g asparagine + 0.5 g glucose (± 0.5 g chitosan, ± 0.5 g HMF) in 1 % acetic acid or in distilled water, pH set to 5.8 with 1 N NaOH then 6.0 with 0.001 N NaOH, topped up to 100 mL, heated at 180 °C for 10, 20 and 30 min, cooled in tap water (sec. 2.2, verbatim in the bundles) |
| What is measured | reducing sugar, asparagine, pH, browning, acrylamide (HPLC-UV, calibration 0–3125 ppb), HMF (HPLC, calibration 0.48–750 ppm), kinematic viscosity |
| Bundles | `mp_holdout_glucose_asparagine_180C_10min_Chang2021`, `..._30min_Chang2021` (acetate arm), `..._30min_water_Chang2021` (water arm) |

## 1. What this paper is and is not, for this model

A **time-course model study** (three times, one temperature) whose acrylamide and HMF levels are in
figures (Figures 3 and 4), with a handful of them repeated in the prose. The bundles score only the
prose-printed numbers and label them "TEXT-QUOTED, FIGURE-SOURCED"; the print confirms each quote
word for word. The heating apparatus is **not named** (the 2022 companion names a dry-bath
incubator); the vessel, its closure and the atmosphere stay unstated.

## 2. Numbers printed in the prose (verbatim)

| quantity | value | where |
|---|---|---|
| acrylamide, 0.5 % glucose + 0.5 % asparagine in 1 % acetic acid at pH 6, 10/20/30 min | **28, 912, 1459 ppb** | sec. 3.3, "(Figure 3B)" |
| acrylamide, the same in deionized water at pH 6, 30 min | **832 ppb** | sec. 3.3, "(Figure 3A,B)" |
| HMF, asparagine + glucose arm, 30 min | **7 ppm** | Discussion, "(Figure 4A)"; the 30 ppm in the same sentence is the glucose-only arm |

## 3. The bundles, checked against the print (wave B36)

| bundle | value | print | verdict |
|---|---|---|---|
| 10 min, acetate | acrylamide 28 ppb | 28 ppb | matches |
| 30 min, acetate | acrylamide 1459 ppb | 1459 ppb | matches |
| 30 min, water | acrylamide 832 ppb; HMF 7000 ppb | 832 ppb; 7 ppm | matches |

The charge (27.75 mM glucose, 33.3 mM asparagine as the monohydrate) is the bundles' own
arithmetic from 0.5 g in 100 mL and is unchanged. The 1 % acetic acid's w/v-versus-v/v ambiguity the
buffer note reports is real in the print too: the paper says "1% acetic acid" and nothing more.
