# Wong & Bernhard 1988 — EXTRACTION
DOI: 10.1021/jf00079a032 (printed: none; Crossref 2026-11-01 resolves to "Effect of nitrogen source on pyrazine formation", J. Agric. Food Chem. 36, 123-129, 1988, authors Wong / Bernhard)
### *"Effect of Nitrogen Source on Pyrazine Formation"*
### D-glucose 1 M + NaOH 1 M + nitrogen source 8 M (ammonium hydroxide, ammonium formate, ammonium acetate, glycine, monosodium glutamate) in water, 120 °C solution temperature, 4 h; end pH ~3 units below the start in every pot. Table II: relative weight % of 32 identified compounds (19 pyrazines), 4 replicates, mean ± SD, FID response factors measured (Table I), capillary columns. **All four target pyrazines resolved, with SDs: R per N source is computable with an error.**

**Read 2026-11-01 from `data/articles/wong1988.pdf`** (7 pp.; the first page carries the end of the preceding JAFC
article on shrimp ketones, pp. 123; Wong & Bernhard is JAFC 36(1):123-129). Acrobat Capture scan with OCR; Tables I-II
from `pdftotext -layout`, cell values spot-checked against the layout (SDs in parentheses under each mean). No data
figure (Fig. 1 is a chromatogram, Figs. 2-3 schemes). Every derived number was computed in the `maillard` env.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed · **[F]** read off a figure · **[Q]** qualitative
claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Authors | Jennifer M. Wong, Richard A. Bernhard (Food Science and Technology, UC Davis) |
| Venue | J. Agric. Food Chem. 36 (1), 123-129, 1988; received 3 Mar 1987, revised 27 Aug 1987, accepted 29 Sep 1987 [M] |
| Laboratory | UC Davis (Bernhard): `ucdavis_shibamoto`. Same lab as shibamoto1976 / 1977, bemisyoung1993 |

## 1. What this paper answers for the model

| need | what Wong 1988 gives | tag |
|---|---|---|
| PZ / MPZ / 2,5-DMP / 2,6-DMP per pot | **yes, five N sources, mean ± SD of 4 replicates** (Table II, relative weight % of identified compounds) | [M] |
| R (molar; factor 0.978) ± 1 SD (dex, quadrature) | **NH4OH 14.6 (± 0.12 dex); NH4 formate 556 (± 0.22); NH4 acetate 349 (± 0.15); glycine 5.4 (± 0.31); glutamate 0.28 (± 0.24)** | [Z] |
| response factors | FID relative RFs measured for every compound (Table I; PZ 1.000, MPZ 1.164, 2,5-DMP 1.177, 2,6-DMP 1.271) by McNair & Bonelli [C]; whether the work-up recovery was checked is not stated here (the method is "described previously", Shibamoto & Bernhard 1978 [C], not on disk) | [M] / [!] |
| 2,5- vs 2,6-DMP | resolved on DB-Wax+ and DB-1701 capillaries; identities by two Kovats indices + MS | [M] |
| pH dependence | **start pH not printed** [!]; every pot dropped ~3 units by the end [M]. NH4OH + NaOH is the most basic start; ammonium salts with 0.1 mol NaOH against 0.8 mol NH4+ start near pH ~8-9 [Z] | [M] / [Z] |
| temperature | one (120 °C, 4 h) | — |
| nitrogen source | **the paper's factor**: five sources at equal molar N | [M] |
| separates P1x / P1y from P2b? | **No** (no labels; no supply perturbation at fixed chemistry) | — |

**Answers.**

1. **R depends on the N source over 3.3 dex at fixed glucose, base, T and time** [Z]: ammonium formate / acetate
   350-560 (MPZ ~45 % of pyrazines), NH4OH 14.6, glycine 5.4, glutamate 0.28. The two ammonium salts (lower start pH)
   agree with UC Davis's glucose + NH3 values of 1976 / 1977 (310-700); NH4OH in this paper does not (disagreement 1).
2. **Glutamate R 0.28 ± 0.24 dex** [Z]: a third low-R row (with the glycine rows of bemisyoung1993 and CF-087) that
   every pairing rival forbids unless h < 0.07 or an extra 2,5-DMP / PZ route exists. Here 2,5-DMP is 18 % of identified
   compounds and acetol 21 % [M]; the authors propose acetaldehyde from glutamate degradation feeding C3 units [P].
3. **Glycine pot: 2,3,5-trimethylpyrazine 28.6 % and tetramethyl 7.1 %** of identified compounds vs MPZ 0.78 % [M]:
   aldehyde (HCHO) addition dominates; the authors attribute extra 2,6-DMP and TrMP to formaldehyde [P].

## 2. Conditions [M]

| item | value |
|---|---|
| Reactants | D-glucose 0.1 mol + NaOH 0.1 mol + N source 0.8 mol in 100 mL deionised water (1 M / 1 M / 8 M) |
| N sources | ammonium hydroxide, ammonium formate, ammonium acetate, glycine, monosodium glutamate |
| pH | measured before and after; values not printed [!]; "dropped approximately 3 pH units for all samples" [V] |
| T / time | 120 °C reaction solution temperature, 4 h ("selected to maximize yield", citing shibamoto1976) |
| Vessel | not described here; method of Shibamoto & Bernhard 1978 [C] (continuous extraction) |
| GC | HP 5711A FID, 30 m DB-Wax+; Varian Vista 6000 FID, 60 m DB-1701; GC-MS VG ZAB-HS-2F, 70 eV |
| Yields | CH2Cl2 extract / glucose: glycine 14.99 %, glutamate 15.64 %, NH4OH 10.38 %, formate 8.65 %, acetate 7.66 % [M] |
| Replicates | 4; mean ± SD |

## 3. Quantitative results (Table II, relative weight % of identified compounds; SD in parentheses) [M]

| compound | NH4OH | NH4 formate | NH4 acetate | glycine | glutamate |
|---|---:|---:|---:|---:|---:|
| pyrazine | 1.137 (0.198) | 0.161 (0.044) | 0.155 (0.024) | 0.016 (0.007) | 0.593 (0.228) |
| 2-methylpyrazine | 14.789 (1.457) | 17.101 (3.638) | 14.392 (2.249) | 0.784 (0.219) | 1.753 (0.328) |
| 2,5-dimethylpyrazine | 12.923 (0.974) | 3.194 (0.261) | 3.749 (0.241) | 6.922 (0.741) | 18.035 (2.064) |
| 2,6-dimethylpyrazine | 17.367 (1.189) | 13.262 (0.759) | 13.964 (1.615) | 3.053 (0.277) | 0.421 (0.195) |
| trimethylpyrazine | 15.055 (1.073) | 1.119 (0.168) | 1.469 (0.108) | 28.591 (0.770) | 2.744 (0.801) |
| tetramethylpyrazine | 3.494 (0.316) | nd | 0.043 (0.006) | 7.137 (0.751) | 0.103 (0.018) |
| pyrrole | 0.932 (0.116) | 0.070 (0.027) | 0.243 (0.027) | 1.075 (0.120) | 0.078 (0.031) |
| acetamide | 5.422 | 7.207 | 44.629 | 0.839 | 0.613 |
| 4-methylimidazole | nd | 38.659 | 7.461 | 2.523 | 0.941 |

Derived [Z] (molar, factor 0.978): R(2,5) = 14.6 / 556 / 349 / 5.4 / 0.28; R(2,5 + 2,6) = 6.2 / 108 / 74 / 3.8 / 0.27;
2,6/2,5 = 1.34 / 4.15 / 3.72 / 0.44 / 0.02. NH4OH pot as a share of the 19 pyrazines (sum 78.4 %): PZ 1.45, MPZ 18.9,
2,5-DMP 16.5, 2,6-DMP 22.2 % [Z].

## 4. Mechanism evidence and its kind

- **Within-study factor (kind 6), one lab:** N source at fixed glucose / base / T / time.
- No labels, no fed intermediates. Glutamate -> 4-oxobutanoic acid -> acetaldehyde + 2-hydroxyacetate (Fig. 3) and
  formaldehyde-derived C2 / C3 fragments in the glycine pot are **proposed** [P].

## 5. Bearing on open questions

| item | effect | why |
|---|---|---|
| EX-052 P0 / P1 | CONTRADICTS P0 and P1 in the ammonium pots (R 14.6-556) | calibrated, resolved isomers, SDs |
| P1x / P1y vs P2b | LEAVES UNSEPARATED | no labels or supply perturbation |
| CF-087 / rivals' floor R >= 4h, h >= 0.1 | CONTRADICTS the floor in a third system (glutamate 0.28) | needs extra 2,5-DMP / PZ routes or strong segregation |
| E4 (pH law) | LEAVES UNCHANGED | start pH not printed |
| evidence map 2.7 (aldehyde addition) | SUPPORTS direction | glycine pot TrMP 28.6 %, TeMP 7.1 % |

## 6. Disagreements found

1. **Same lab, glucose + NH4OH: R 14.6 (± 0.12 dex) here vs 310-700 in shibamoto1976 / 1977** (1 M glucose + 8 M NH3;
   1976 rows at 100 °C 2-48 h, 120-160 °C closed, ± 0.1 M NaOH). 1.3-1.7 dex. The largest change is 2,5-DMP relative to
   PZ (11.4 here vs ~1.3 there) and MPZ's share of pyrazines (19 vs ~80 %). Candidate causes: (a) NaOH 1 M here vs
   0-0.1 M there (1976 exp 43-44 with 0.1 M NaOH still give R 430-530, so a 10x larger dose would have to act
   non-linearly); (b) work-up: continuous extraction (Shibamoto & Bernhard 1978) vs four batch extractions; the authors
   blame "the method of extraction and separation capabilities" [V] for the differences with 1977;
   (c) packed (1976 / 1977) vs capillary columns, i.e. 2,5-DMP co-eluting with or misassigned against 2,6-DMP /
   2-ethylpyrazine (Kovats 1320 / 1326 / 1330 on Wax); (d) 120 vs 100 °C (1976's 120 °C row duplicates its 100 °C row,
   so T is untested there, see shibamoto1976 flag).
2. **Ammonium salts vs NH4OH in this paper: R 350-560 vs 14.6** — opposite in direction to the glycine pH trend of
   bemisyoung1993 (R up with pH). Candidate causes: start pH (salts lower), acetamide / formate side chemistry (acetate
   sequesters N as acetamide, 44.6 %).
3. **Glycine MPZ** 0.78 % of identified compounds here vs koehler1970 ND (glycine-glucose, glycol, no base).

## 7. Ledger-ready measurements

1. Table II, PZ / MPZ / 2,5-DMP / 2,6-DMP, five N sources, mean ± SD (n = 4) — the only calibrated, isomer-resolved
   N-source ladder of R on disk; key as `relative_weight_percent_of_identified_rf_corrected`.
2. Table I response factors (FID, relative to pyrazine) — usable to convert other UC Davis area tables.
3. Glycine column TrMP / TeMP — aldehyde-addition products (direction).
