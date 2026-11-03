# Matoba, Kuchiba, Kimura & Hasegawa 1988 — EXTRACTION (file stem `matoba1988`)
DOI: 10.1111/j.1365-2621.1988.tb13551.x (not printed on the pages; Crossref 2026-11-03: "Thermal Degradation of
Flavor Enhancers, Inosine 5′-monophosphate, and Guanosine 5′-monophosphate in Aqueous Solution", J. Food Sci.
**53**(4):**1156-1159**, authors MATOBA, KUCHIBA, KIMURA, HASEGAWA, issued 1988-07. **Matches the PDF and the brief.**)
[!] **The PDF is incomplete**: p. 1159 ends "—Continued on page 1170"; p. 1170 (end of the mechanism discussion,
any conclusion, the reference list) is **not in the file** (4 pp., 1156-1159). Nothing tabular is lost (Tables 1-3
and Figs. 1-6 are all on 1156-1159), but the cited references cannot be checked.
**Companion paper, same group:** `kuchiba1990_extraction.md` (salt effects, same buffers, same HPLC).

### IMP or GMP 0.4 mM in 0.01 M buffer (acetate pH 4.0, MOPS pH 7.0, CHES pH 9.0), sealed Pyrex tubes, 80-100 °C, 0-7.5 h (kinetics) and 0-200 h (pH 4, 100 °C, products). First-order k ± for IMP and GMP at 3 pH × 2-3 T, printed half-lives and Ea (Tables 1-2); Pi-release k for ribose 5-phosphate, O-phosphorylethanolamine and glycerophosphate (Table 3); time courses of nucleotide, nucleoside and purine base by HPLC-UV (Figs. 2-4). **No ribose, no R5P measured from the nucleotides; products only at pH 4; no T above 100 °C; no pH after heating; no replicates stated.**

**Read 2026-11-03 from `data/articles/matoba1988.pdf`** (4 pp., scan with text layer). Prose from `pdftotext
-layout`; **all four pages rendered at 200 dpi** (scratchpad `ingest/matoba1988_work/p-1..4.png`) and Tables 1-3,
Figs. 2-5 read from the renders; **Fig. 2 re-read from a 300 dpi crop** (`fig2b-2.png`), read error ±8 µM on the
0-400 µM axis and ±2 h on the 0-200 h axis; Figs. 3-4 read from the 200 dpi render, ±8 µM, ±0.2 h. Arrhenius
re-fits of Tables 1-3 by `awk` on `matoba1988_work/k.tsv` [Z].

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Dept. of Food Science & Nutrition, Faculty of Home Economics, **Nara Women's University**, Nara
(Matoba, Kuchiba, Hasegawa); Kimura at the Dept. of Chemistry, Faculty of Science, same university. No
corresponding author marked; **Teruyoshi Matoba** is the corresponding author of the companion Kuchiba 1990 and
Hasegawa is the senior author on both. **Not in `data/keys/laboratories.yml`: NEW**, proposed id
**`nara_womens_matoba`**, lab string "Nara Women's Univ. (Matoba / Hasegawa)" (same proposal in
`kuchiba1990_extraction.md`). No other dossier from this group; `zhang1991_extraction.md` (Zhang & Ho 1991, IMP +
Cys at 180 °C) cites this paper for IMP hydrolysis but is Rutgers, unrelated.

## 1. What this paper answers for the model

| model need | what Matoba 1988 gives | tag |
|---|---|---|
| IMP / GMP loss rate vs pH and T (aqueous) | Tables 1-2: k_obsd (10⁻² h⁻¹) ±, pH 4 / 7 / 9, 80-100 °C (2-3 T per pH) | [M] |
| Ea | per pH, 1.1-1.6 × 10⁵ J/mol; "about 1.3 × 10⁵ J·mol⁻¹ (30 kcal·mol⁻¹)" overall | [M] |
| pre-exponential A | **not printed**; my ln A per pH (sec. 4b) | [Z] |
| half-lives | printed, t½ = 0.693/k (Tables 1-2, abstract) | [M] |
| product split nucleoside + Pi vs base + R5P | two pathways drawn (Fig. 6); split not printed; my estimate from Figs. 3-4 (~18 % IMP, ~5 % GMP via base + R5P at pH 4) | [P] / [Z] |
| free ribose / R5P yield | **not measured** (inferred from free base by mass balance) | [!] |
| nucleoside N-glycosidic cleavage k | text: 1.8 × 10⁻² h⁻¹ inosine, 0.6 × 10⁻² h⁻¹ guanosine (pH 4, 100 °C, from Figs. 3-4) | [M] |
| R5P dephosphorylation k | Table 3: 8.3 / 3.2 / 1.1 × 10⁻² h⁻¹ at 100 / 90 / 80 °C, pH 4; Ea 1.2 × 10⁵ | [M] |
| T > 100 °C (retort, extrusion, roasting) | **none** | [!] |
| pH during / after heating; ionic strength | **not reported**; 0.01 M buffers, no salt | [!] |
| matrix (meat, yeast extract, protein) | **none**; dilute buffer only | [!] |

**Answers.**
1. **Rate constants vs pH and T: yes, clean first-order k** (confidence ~90 % that the transcription in sec. 4a is
   exact; ~80 % that the k are good descriptions of the 0-7.5 h loss). Units 10⁻² h⁻¹. At 100 °C: IMP 8.0 / 5.3 / 1.5,
   GMP 11.0 / 8.5 / 1.8 (pH 4 / 7 / 9). Acid faster; GMP faster than IMP at every pH [M]. Half-lives and k agree
   (t½ = 0.693/k) in every cell to rounding except GMP pH 9 / 90 °C (0.693/0.006 = 115.5 vs 108.3 printed: the
   printed k 0.6 is rounded from ~0.64) [Z].
2. **Ea / A: Ea printed and reproducible; A not printed** (confidence ~90 %). My three-point fits give IMP 114 /
   128 / 174 (two points) kJ/mol at pH 4 / 7 / 9 and GMP 117 / 112 / 124 — printed 1.1 / 1.3 / 1.6 and 1.2 / 1.1 /
   1.2 × 10⁵ J/mol [Z]. IMP pH 9 (174) rounds to 1.7, printed 1.6; it rests on two points 5 K apart, so treat as
   poorly determined (± ~30 kJ/mol). "One-third per 10 °C" (abstract) is a Q10 ≈ 2.7-3.0 over 80-100 °C [Z].
3. **Half-lives: yes** [M], at 100 °C IMP 8.7 / 13.1 / 46.2 h, GMP 6.4 / 8.2 / 38.5 h (pH 4 / 7 / 9).
4. **Products** (confidence ~85 %): nucleoside + phosphoric acid is the main route (pathway 1, "Fast" ester
   hydrolysis, then "Slow" N-glycosidic cleavage to base + ribose); base + ribose 5-phosphate directly from the
   nucleotide is the minor route (pathway 2, "Slow" then "Fast" R5P → ribose + Pi) [P, Fig. 6]. Hypoxanthine,
   inosine, guanine, guanosine are measured; **ribose, R5P and (for the nucleotides) Pi are not shown** [!].
   Mass balance nucleotide + nucleoside + base = 100 % to 50 h, ~95 % after 100 h, with a trace unidentified HPLC
   peak [Q]. My estimate of the pathway-2 share at pH 4 / 100 °C from Fig. 3a / 4a at 7.5 h (base formed minus the
   base expected from nucleoside cleavage at the text's nucleoside k): **IMP ~18 % (15-22 %), GMP ~5 % (3-8 %)** of
   nucleotide lost [Z] (~60 % confidence in the ranges; pH 4 only).
5. **Buffer and ionic strength** [M]: 0.01 M sodium acetate (pH 4.0), 0.01 M MOPS (pH 7.0), 0.01 M CHES (pH 9.0);
   no salt; ionic strength not stated, ≤ ~0.01 M [Z]. Counter-ion of IMP / GMP (likely disodium salts) not stated.
6. **What this means for IMP as a ribose source in cooking** [Z, ~70 % confidence in the order of magnitude]:
   interpolating ln k linearly in pH gives k(IMP, pH 6, 100 °C) ≈ 0.061 h⁻¹ (t½ ≈ 11 h), so 30 min of boiling loses
   ~3 % of IMP (GMP ~4.5 %), mostly to inosine, which keeps its ribose bound; free base (= ribose + R5P released) at
   0.5 h is ~1 % of initial (Fig. 3a, ±1 pp). Extrapolating with the pH 7 Ea (128 kJ/mol) gives k(120 °C) ≈ 0.43 h⁻¹
   (t½ ≈ 1.6 h) and k(140 °C) ≈ 2.9 h⁻¹ — 20-40 K beyond the data, low confidence (~40 % within ×2). **IMP is a
   slow ribose (R5P) source at 100 °C; it matters only for long or pressurised heating.**

## 2. Conditions

- **Reactants** [M]: IMP, GMP, R5P, O-phosphorylethanolamine, β-glycerophosphoric acid (Nakarai Chemicals, Kyoto);
  inosine, hypoxanthine, guanosine, guanine standards (Wako). Purity and salt form not stated.
- **Concentration** [M]: nucleotides **0.4 × 10⁻³ mol dm⁻³** (Tables 1-2 footnotes; Fig. 5 ln[µM] starts at ~6.0 =
  400 µM [F]); esters 2-3 × 10⁻³ mol dm⁻³ (Table 3).
- **Buffer / pH** [M]: as answer 5; pH before heating only; **no pH after heating** [!]. Pi released (≤ 0.4 mM)
  is small against 10 mM buffer [Z].
- **Vessel / heating** [M]: sealed Pyrex tube, incubated at 80-100 °C (heating medium not stated); quench in ice
  bath. Headspace / O2: not stated (air presumably). Heat-up not reported.
- **Time** [M]: kinetics 0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5 h; products at pH 4 / 100 °C to 200 h (Fig. 2).
- **Replicates**: **not stated**; the ± on k undefined (fit SE or replicate SD) [!].

## 3. Analytical methods

- **HPLC-UV** [M]: Shimadzu LC-6A, 254 nm, Cosmosil 5C18-P (4 × 250 mm), methanol-0.02 M KH2PO4 (15:85), 0.5
  mL/min, 25-30 °C; peak areas by C-R3A integrator. Retention: nucleotide 6.5, base 8.5, nucleoside 10.5 min.
  Quantification against purchased standards presumably (calibration not described) [!]. No LOD.
- **TLC** [M]: Kieselgel 60 F254, ethanol-1 M ammonium acetate 2:1, UV detection (Fig. 1, identification only).
- **Pi** [M]: Fiske & Subbarow 1925 molybdate method — used for the esters (Table 3); not shown for IMP / GMP.
- **Kinetics** [M]: slope of ln[nucleotide] vs t (Fig. 5); Arrhenius ln k = −Ea/RT + ln C. Linear "until both
  IMP and GMP were completely degraded" [Q] (but see flag 2).

## 4. Data

### 4a. Tables 1-2: first-order k_obsd (×10⁻² h⁻¹) ± (undefined), t½ (h), Ea (×10⁵ J mol⁻¹) [M]

| pH | T (°C) | IMP k | IMP t½ | GMP k | GMP t½ |
|---|---|---|---|---|---|
| 4 | 100 | 8.0 ± 0.12 | 8.7 | 11.0 ± 0.11 | 6.4 |
| 4 | 90 | 2.7 ± 0.12 | 25.7 | 3.7 ± 0.03 | 18.7 |
| 4 | 80 | 1.0 ± 0.05 | 69.3 | 1.3 ± 0.04 | 53.3 |
| 7 | 100 | 5.3 ± 0.002 | 13.1 | 8.5 ± 0.02 | 8.2 |
| 7 | 95 | 2.9 ± 0.06 | 23.9 | — | — |
| 7 | 90 | 1.7 ± 0.09 | 40.8 | 3.0 ± 0.01 | 23.1 |
| 7 | 80 | — | — | 1.1 ± 0.06 | 63.0 |
| 9 | 100 | 1.5 ± 0.08 | 46.2 | 1.8 ± 0.09 | 38.5 |
| 9 | 95 | 0.7 ± 0.05 | 99.0 | 1.0 ± 0.07 | 69.3 |
| 9 | 90 | — | — | 0.6 ± 0.05 | 108.3 |

Ea (×10⁵ J mol⁻¹): IMP 1.1 (pH 4), 1.3 (pH 7), 1.6 (pH 9); GMP 1.2, 1.1, 1.2. [!] The GMP table's printed title
says "inosine 5′-monophosphate (IMP)" (it is Table 2; footnote a says GMP solution) — a typesetting error.

### 4b. My Arrhenius re-fit (unweighted ln k vs 1/T, k in h⁻¹) [Z]

| series | n | Ea (kJ/mol) | ln A (A in h⁻¹) | k(100 °C) fit | k(120 °C) extrap. |
|---|---|---|---|---|---|
| IMP pH 4 | 3 | 114 | 34.14 | 0.078 | 0.50 |
| IMP pH 7 | 3 | 128 | 38.33 | 0.052 | 0.43 |
| IMP pH 9 | 2 | 174 | 51.92 | 0.015 | 0.26 |
| GMP pH 4 | 3 | 117 | 35.47 | 0.108 | 0.74 |
| GMP pH 7 | 3 | 112 | 33.61 | 0.084 | 0.52 |
| GMP pH 9 | 3 | 124 | 35.85 | 0.018 | 0.14 |
| R5P pH 4 (Table 3) | 3 | 111 | 33.22 | 0.084 | 0.52 |

### 4c. Table 3: Pi release from phosphoric esters, k_obsd (×10⁻² h⁻¹) ±, Ea (×10⁵ J mol⁻¹) [M]

| ester | pH | 100 °C | 95 °C | 90 °C | 80 °C | Ea |
|---|---|---|---|---|---|---|
| A ribose 5′-phosphate | 4 | 8.3 ± 0.12 | — | 3.2 ± 0.04 | 1.1 ± 0.02 | 1.2 |
| B O-phosphorylethanolamine | 4 | 10.4 ± 0.06 | — | 3.6 ± 0.003 | 1.1 ± 0.02 | 1.3 |
| B | 7 | 3.3 ± 0.03 | 1.8 ± 0.01 | 1.1 ± 0.02 | — | 1.3 |
| C glycerophosphoric acid | 4 | 9.6 ± 0.02 | — | 3.0 ± 0.07 | 0.9 ± 0.04 | 1.3 |
| C | 7 | 6.8 ± 0.04 | — | 2.0 ± 0.004 | 0.5 ± 0.08 | 1.6 |
| C | 9 | 1.5 ± 0.01 | 0.9 ± 0.01 | — | — | 1.2 |

R5P at pH 7 / 9 **not measured**. Text: nucleoside N-glycosidic cleavage at pH 4 (100 °C, from Figs. 3-4) **1.8 ×
10⁻² h⁻¹ inosine, 0.6 × 10⁻² h⁻¹ guanosine** [M].

### 4d. Fig. 2 (100 °C, pH 4, 0-200 h), µM, read ±8 µM, times ±2 h [F]

| t (h) | IMP | inosine | Hx | | t (h) | GMP | guanosine | guanine |
|---|---|---|---|---|---|---|---|---|
| 0 | 400 | 0 | 0 | | 0 | 373 | 0 | 0 |
| 4 | 314 | 66 | 19 | | 3.5 | 238 | 105 | 10 |
| 9 | 187 | 155 | 49 | | 8.5 | 124 | 224 | 24 |
| 24 | 24 | 213 | 133 | | 24 | 37 | 297 | — |
| 32 | 11 | 175 | 191 | | 32 | 10 | 285 | 84 |
| 53 | ~0 | 112 | 266 | | 52 | 6 | 231 | 139 |
| 72 | ~0 | 76 | 297 | | 71 | ~2 | 184 | 183 |
| 80 | ~0 | 44 | 326 | | 79 | ~2 | 163 | 192 |
| 96 | ~0 | 23 | 343 | | 95 | ~2 | 131 | 219 |
| 119 | ~0 | 8 | 356 | | 142 | ~1 | 92 | 280 |
| 143 | ~0 | 6 | 361 | | 169 | ~0 | 58 | 310 |
| 193 | ~0 | 6 | 360 | | 190 | ~0 | 51 | 320 |

Inosine peaks ~215 µM at ~25 h; Hx plateaus ~360 µM (~90 % of initial) after ~120 h; guanosine peaks ~300 µM at
~24 h. [Z] Mass balance 92-95 % from 24 h on (both panels). Tail slopes after the nucleotide is gone:
inosine ln(175/23)/64 h ≈ 0.032 h⁻¹, guanosine ln(285/51)/158 h ≈ 0.011 h⁻¹ — **~1.8 × the text's 7.5 h values**
(flag 3).

### 4e. Figs. 3-4 (100 °C, pH 4, 0-7.5 h), µM at 7.5 h, read ±8 µM [F]

| start | nucleotide left | nucleoside | base |
|---|---|---|---|
| IMP 400 (Fig. 3a) | 222 | 128 | 40 |
| inosine 400 (Fig. 3b) | — | 350 | 57 |
| GMP 400 (Fig. 4a) | 180 | 185 | 15 |
| guanosine 400 (Fig. 4b) | — | 378 | 18 |

[Z] ln(400/222)/7.5 = 0.079 h⁻¹ and ln(400/180)/7.5 = 0.106 h⁻¹ reproduce Tables 1-2; ln(400/350)/7.5 = 0.018 and
ln(400/378)/7.5 = 0.0075 reproduce the text's nucleoside k. Text: Hx from IMP at 7.5 h is "about 60 %" of that from
inosine [Q] (my read 40/57 = 70 %). Fig. 1 (TLC, 0-3 h) and Fig. 5 (ln plot) are qualitative / duplicate.

## 5. Flags and internal consistency

1. [!] **p. 1170 missing** from the PDF (continuation + references).
2. [!] **Fig. 2 IMP falls faster than Table 1 after 7.5 h**: 24 µM at ~24 h vs 400·e^(−0.080·24) = 59 µM predicted
   (beyond the ±8 µM read error); GMP 37 vs ~28 predicted (fine). The "first order until completely degraded" claim
   holds for GMP, not clearly for IMP.
3. [!] Nucleoside N-glycosidic k from the long Fig. 2 tails (~0.032 / 0.011 h⁻¹) is ~1.8 × the 7.5 h text values
   (0.018 / 0.006) [Z]; the short-run values rest on 5-12 % loss and are the less precise.
4. [!] Table 2 title misprinted (IMP for GMP). GMP pH 9 / 90 °C t½ 108.3 implies k = 0.64, printed 0.6.
5. [!] IMP pH 9 Ea from two temperatures 5 K apart (1.6, my fit 1.74 × 10⁵); R5P Ea 1.2 printed vs my 1.11.
6. [!] Fig. 2b GMP at t = 0 reads ~373 µM, not 400.
7. [!] The ± on k is undefined; IMP pH 7 / 100 °C ± 0.002 (0.04 % RSD) is implausibly small. No replicates stated.
8. [!] Fig. 6 labels R5P → ribose + Pi "Fast", yet Table 3 R5P k (8.3 × 10⁻²) ≈ IMP k (8.0 × 10⁻²) at pH 4 / 100 °C:
   R5P is about as long-lived as IMP, so an R5P pool can accumulate. R5P's own Maillard / degradation chemistry
   (not Pi-releasing) is invisible to the Pi assay.
9. [!] Text: "Ea … about 1.3 J·mol⁻¹" for the esters (p. 1158) is a missing ×10⁵ (Table 3 has it).

## 6. Existing records in the repo that cite this paper

Grep 2026-11-03 (`grep -rIn -iE 'matoba|kuchiba|tb13551|nucleotide_pathway_priors|00092a'` over `data docs src
scripts tasks`, plus `results`, `tests`). Records: `computational_priors.json` → `nucleotide_pathway_priors` →
`matoba_1988_nucleotide_hydrolysis_v1` and `matoba_1988_nucleotide_hydrolysis`; `benchmark_intake_registry.json`
→ `matoba_1988_nucleotide_hydrolysis` (l. 524); `data/keys/papers.yml` l. 7386; `data/research_corpus/raw/
Gemini_deep_research_04.md` l. 129-148; `tasks/data_restructure_plan.md` l. 1168 (D8); `docs/guides/
LITERATURE_SEARCH_2026-10_ROUND3.md` l. 37 / 350; `results/validation/key_value_trace_report.md` l. 129-130;
`citation_verification_ledger.md` l. 606. **No "Kuchiba" in `computational_priors.json`** (grep: none).

### 6a. `computational_priors.json` → `matoba_1988_nucleotide_hydrolysis_v1`

| field | repo value | paper | verdict |
|---|---|---|---|
| source | "Matoba, Terao & Fujimaki 1988 … plus Mouritsen 2024 …, Yamaguchi & Ninomiya 2000 …, Nakamura 1988 …" | Matoba, Kuchiba, Kimura & Hasegawa 1988 | **wrong co-authors** (no Terao, no Fujimaki); the other three sources are not this paper |
| doi / doi_propagation.crossref_verified | tb13551.x; MATOBA, KUCHIBA, KIMURA, HASEGAWA, 1988, JFS | same | **matches** (and contradicts the record's own source string) |
| yeast_extract_imp_share / gmp_share | 0.6 / 0.4 | no yeast extract | **not in paper** |
| umami_weight / ribose_weight | 0.6 / 0.4 | — | **not in paper** (repo weights) |
| imp / gmp_fold_reduction_at_3mm | 45.2 / 29.8 | no sensory data | **not in paper** (attributed to Mouritsen) |
| imp / gmp_relative_umami_unit, umami_synergy_constant | 1.0 / 2.3, 1218 | — | **not in paper** (Yamaguchi EUC) |
| imp_half_life 100 °C pH 4.0 / 7.0 / 9.0 | 8.7 / 13.1 / 46.2 | 8.7 / 13.1 / 46.2 | **match** |
| imp_half_life 100 °C pH **5.0** | 12.0 | no pH 5 run | **not in paper** |
| gmp_half_life 100 °C pH 4.0 / 7.0 / 9.0 | 6.4 / 8.2 / 38.5 | 6.4 / 8.2 / 38.5 | **match** |
| gmp_half_life 100 °C pH **5.0** | 9.0 | no pH 5 run | **not in paper**, and **contradicts its trend** (9.0 > 8.2 at pH 7; the paper has t½ rising monotonically with pH) |
| imp_half_life pH 7 at 100 / 110 / 130 / 150 °C | 13.1 / 4.4 / 0.5 / 0.0583 | 100 °C only = 13.1; data stop at 100 °C | 100 **matches**; 110 / 130 / 150 **not in paper** — "one-third per 10 °C" extrapolated (13.1/3 = 4.37, /27 = 0.49, /243 = 0.054); my Ea-128 Arrhenius gives 4.5 / 0.61 / **0.10** h |
| gmp_half_life pH 7 at 100 / 110 / 130 / 150 °C | 8.2 / 2.7 / 0.3 / 0.0367 | 100 °C = 8.2 | 100 **matches**; rest **not in paper** (Arrhenius with GMP pH 7 Ea 112 gives 150 °C ≈ 0.12 h, 3 × the repo) |
| default_extrusion_residence_minutes, ribose_release_efficiency_multiplier, max_ribose_delivery_fraction, ribose_addition_threshold, ribose_addition_ratio_scale | 2.0, 2.3, 0.4, 0.08, 0.35 | — | **not in paper** |
| reference_ribose_release_fraction_100c_30min | 0.06 | total IMP loss in 30 min at 100 °C: 3.9 % (pH 4), 2.6 % (pH 7); free base (= ribose + R5P released) ~1 % (Fig. 3a) [Z] | **not in paper; exceeds even the total IMP loss** (~5-10 × the ribose + R5P release the data imply) |
| notes ("extrusion-like severity") | — | 80-100 °C, dilute buffer | **not in paper** |
| provenance_tier / confidence_tier | literature_derived_transfer / medium | — | tier fine for the half-lives; the record mixes in five other sources' numbers and invented parameters |

### 6b. `computational_priors.json` → `matoba_1988_nucleotide_hydrolysis` (twin, no DOI)

| field | repo value | paper | verdict |
|---|---|---|---|
| source | "Matoba et al. (1988)" | — | **matches** (no DOI on this twin) |
| confidence_tier | high | — | **not supported** by the numbers below |
| reference_conditions | T 80 °C, pH 6.0 | 80 °C exists (pH 4 IMP, pH 4 / 7 GMP); **no pH 6** | **pH not in paper** |
| IMP_hydrolysis_k_1_per_min | 0.0035 (= 0.21 h⁻¹) | IMP 80 °C pH 4: 1.0 × 10⁻² h⁻¹ = **1.7 × 10⁻⁴ min⁻¹**; pH 6 would be slower | **mismatch, ~21 × too fast** |
| IMP_to_inosine_Ea_kj_mol | 76.0 | IMP 110-160 kJ/mol (1.1-1.6 × 10⁵ J/mol); my fits 114-174 | **mismatch** (no 76 anywhere) |
| notes "IMP hydrolysis to inosine" | — | k is total IMP loss (both pathways), not inosine formation | **imprecise** |

### 6c. `benchmark_intake_registry.json` → `matoba_1988_nucleotide_hydrolysis` (l. 524)

| field | repo value | paper | verdict |
|---|---|---|---|
| citation / doi | Matoba et al. (1988) / tb13551.x | same | **match** |
| matrix_family | aqueous_model_system | dilute buffer | **match** |
| process_state_scope | heated_matrix, aqueous_pre_extrusion_model, extrusion_structured | 80-100 °C buffer | **extrusion_structured: not in paper** |
| key_values half-life dicts (by pH, by temp) | as 6a | — | same verdicts as 6a: pH 5.0 and 110-150 °C **not in paper** |
| key_values.replicates | 3 | **not stated** | **not in paper** |
| key_values.assay_mode | hplc_uv | HPLC, UV 254 nm | **match** |
| what_it_supports[0] "survival under heating and extrusion-like severity" | — | heating to 100 °C only | **partly** (no extrusion data) |
| what_it_supports[2] separate IMP and GMP | — | yes | **match** |
| what_it_does_not_support "complete ribose liberation kinetics after inosine cleavage" | — | ribose never measured; nucleoside k at pH 4 only | **correct** |

### 6d. Other mentions

- `Gemini_deep_research_04.md` l. 129-148: the six half-lives **match**; "+10 °C from 100 °C to 110 °C shortens to
  one-third" **misreads the paper** (the one-third is observed over 80-100 °C; nothing above 100 °C) — likely the
  origin of the repo's 110-150 °C values; IMP Ea "34.0 / 30.4 / 28.1 kcal/mol at pH 3 / 4 / 5" is cited to its ref.
  42, **not in Matoba** (pH 4 here ≈ 1.1 × 10⁵ J/mol ≈ 26-27 kcal/mol).
- `key_value_trace_report.md` l. 129 lists 0.0035, 0.0367, 0.0583, 0.3, 0.5, 2.7, 4.4, 76.0, 80 as this record's
  values: all except 80 (°C) are **not in the paper**.
- `papers.yml` l. 7386: `dossier: null` — this file is now the dossier (not edited, per brief).
- `tasks/data_restructure_plan.md` D8 ("IMP hydrolysis kinetics missing: Matoba 1988 PDF") — now on disk.

### 6e. Code that reads these records

`src/kinetic_core/` has no IMP / nucleotide species and no reference to these ids. `src/matrix_prior_registry.py`
loads the whole priors file; its `summarize_family_prior_bundle` / `query_family_prior_entries` would pick up
`nucleotide_and_ribose_support` rows but **have no caller outside the module** (grep). `src/deep_research_runtime_
queue.py` l. 266 harvests every record's `source` string (bookkeeping). `src/literature_runtime.py`, named in the
`doi_propagation` note as the reader, **no longer exists**. **No shipped number depends on these records**
(~90 % confidence).

## 7. Proposed ledger entries (trunk.yml; NOT written)

Trunk is the closest ledger (sugar-side kinetics; IMP / R5P as a ribose donor). None is simulable today (no IMP,
inosine, R5P species).

```yaml
  # ---- matoba1988 ----
  - id: matoba1988_tab1_tab2_imp_gmp_first_order_k_ph_temperature
    dossier: matoba1988_extraction.md
    paper: "Matoba, Kuchiba, Kimura & Hasegawa 1988 J Food Sci 53:1156"
    lab: "Nara Women's Univ. (Matoba / Hasegawa)"
    system: "IMP or GMP 0.4 mM, no amine"
    medium: "0.01 M sodium acetate pH 4.0 / 0.01 M MOPS pH 7.0 / 0.01 M CHES pH 9.0; sealed Pyrex tubes"
    temperature_c: [80, 90, 95, 100]
    ph: {initial: [4.0, 7.0, 9.0], measured_during: false, final: null}
    time_min: [0, 450]
    species: [IMP, GMP]
    kind: rate_constant
    location: "Tables 1-2 (k_obsd x10^-2 h^-1 +/-, t1/2, Ea per pH; dossier sec. 4a-4b)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: "0.01 M sodium acetate at pH 4.0, 0.01M MOPS at pH 7.0"
    flags: [ph_start_only, replicates_not_stated, plus_minus_undefined, table2_title_misprint, no_temperature_above_100c]
    note: "16 k cells (2-3 T per pH); Ea 1.1-1.6e5 J/mol reproduce from the table; A not printed (dossier refit). Repo priors' pH 5 and 110-150 C half-lives are not here."
  - id: matoba1988_tab3_phosphate_ester_hydrolysis_k
    dossier: matoba1988_extraction.md
    paper: "Matoba, Kuchiba, Kimura & Hasegawa 1988 J Food Sci 53:1156"
    lab: "Nara Women's Univ. (Matoba / Hasegawa)"
    system: "ribose 5-phosphate, O-phosphorylethanolamine or glycerophosphoric acid 2-3 mM"
    medium: "0.01 M acetate pH 4 / MOPS pH 7 / CHES pH 9; sealed Pyrex tubes"
    temperature_c: [80, 90, 95, 100]
    ph: {initial: [4.0, 7.0, 9.0], measured_during: false, final: null}
    time_min: [0, 450]
    species: [ribose 5-phosphate, phosphate]
    kind: rate_constant
    location: "Table 3 (Pi release by Fiske-Subbarow; R5P at pH 4 only; dossier sec. 4c)"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: ""
    flags: [ph_start_only, pi_release_only_not_r5p_loss, r5p_ph4_only]
    note: "R5P dephosphorylation k 8.3/3.2/1.1e-2 h^-1 at 100/90/80 C, pH 4; Ea 1.2e5. The Pi assay misses non-hydrolytic R5P loss."
  - id: matoba1988_fig2_fig4_nucleotide_nucleoside_base_time_courses
    dossier: matoba1988_extraction.md
    paper: "Matoba, Kuchiba, Kimura & Hasegawa 1988 J Food Sci 53:1156"
    lab: "Nara Women's Univ. (Matoba / Hasegawa)"
    system: "IMP, GMP, inosine or guanosine 0.4 mM"
    medium: "0.01 M sodium acetate pH 4.0; sealed Pyrex tubes"
    temperature_c: [100]
    ph: {initial: [4.0], measured_during: false, final: null}
    time_min: [0, 12000]
    species: [IMP, inosine, hypoxanthine, GMP, guanosine, guanine]
    kind: time_series
    location: "Fig. 2 (0-200 h) and Figs. 3-4 (0-7.5 h), digitised dossier sec. 4d-4e; nucleoside N-glycosidic k in text p. 1159"
    status: digitised_in_dossier
    quant_class: calibrated
    simulable: false
    condition_quote: ""
    flags: [figure_read_error_8uM, ribose_r5p_not_measured, imp_faster_than_first_order_after_7h, ph4_only]
    note: "Two-pathway split (base + R5P direct ~18 % IMP / ~5 % GMP, dossier estimate) and nucleoside cleavage k 1.8e-2 (inosine) / 0.6e-2 (guanosine) h^-1; long-run tails ~1.8x faster."
```

## 8. Supplementary material

None referenced (1988 J. Food Sci.). The only missing item is **p. 1170** of the same issue (continuation of the
discussion and the reference list, J. Food Sci. 53(4)); nothing quantitative is expected there (~85 %), but it
should be fetched to check the citations (Kuriyama 1965, Fujita 1965, Davidek 1972, Lee 1979, Shaoul & Sporns 1987).
