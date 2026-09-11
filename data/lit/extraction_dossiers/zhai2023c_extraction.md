# Zhai, Xia, Deng, Cui, Hayat, Zhang & Ho 2023 — EXTRACTION (TTCA ± xylose at 100/120/140 °C: 2-methyl-3-furanthiol and 2-furfurylthiol measured directly against time)

**Source on disk:** `data/articles/Zhai2023c.pdf` (1.6 MB; downloaded 2026-09-11). Read 2026-09-11
via `pdftotext -layout` with Table 1 and Fig. 1a verified against 150–400 dpi page renders. Wave B37.

| field | value |
|---|---|
| Title | "Reduced asynchronism between regenerative cysteine and fragments of deoxyosones promoting formation of sulfur-containing compounds through extra-added xylose and elevated temperature during thermal processing of 2-threityl-thiazolidine-4-carboxylic acid" |
| Venue | Food Chemistry 404 (2023) 134420 |
| DOI | 10.1016/j.foodchem.2022.134420 |
| Group | Jiangnan University with Miami University and Rutgers (Chi-Tang Ho) |
| Systems | **TTCA (10 mmol/L) with and without extra xylose (10 mmol/L)**, pH set to 7.0 with NaOH in deionized water, no buffer, pressure-resistant bottles, **100, 120 and 140 °C**, times 20–140 min; plus a cysteine-only control and a fresh xylose–cysteine MRP control (0.0827 mol/L each, 120 °C, 120 min) |
| Quantification | HS-SPME (75 µm CAR/PDMS/DVB), 3 g sample + 2 g saturated NaCl, 60 °C, 20 min; internal standard 1,2-dichlorobenzene, **"the correction factor was 1"** |

## 1. Why this paper matters to the sulfur lane

TTCA is already a chargeable precursor in this model (`PRECURSOR_ALIASES`, wave W6: the
xylose–cysteine thiazolidine, about 94 % of the group's "Cys-Amadori"). This paper runs **the fed
TTCA pot the repository's own species table was built for**, at three temperatures, with a printed
time course, and with **2-methyl-3-furanthiol and 2-furfurylthiol quantified as free thiols** rather
than inferred from their dimers. The thiol-sink question (waves B25 and B27) has been blocked for
want of exactly this.

## 2. Fig. 1a — the time course at 100 °C, µg/L (every cell carries its printed number)

Columns: T = TTCA alone, T-X = TTCA + xylose; the last column O is the fresh xylose–cysteine MRP
(120 °C, 120 min). **No standard deviations are printed in this figure.** Six rows reconcile exactly
with Table 1's 100 °C columns, which fixes Table 1's unstated time as 120 min.

| compound | T-20 | T-X-20 | T-40 | T-X-40 | T-60 | T-X-60 | T-80 | T-X-80 | T-100 | T-X-100 | T-120 | T-X-120 | T-140 | T-X-140 | MRPs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **2-Methyl-3-furanthiol** | 0 | 0.083 | 0.078 | 0.213 | 0.100 | 0.334 | 0.509 | 0.939 | 0.980 | 1.537 | 1.237 | 1.729 | 1.179 | 1.928 | 0.477 |
| **2-Furfurylthiol** | 0.253 | 0.372 | 0.708 | 1.218 | 0.842 | 1.397 | 1.206 | 2.101 | 3.022 | 3.499 | 3.734 | 4.736 | 3.029 | 4.875 | 0.064 |
| 3-Mercapto-2-butanone | 0 | 0 | 0 | 0.028 | 0.079 | 0.129 | 0.198 | 0.349 | 0.239 | 0.386 | 1.303 | 1.673 | 1.529 | 2.018 | 0.728 |
| 3-Mercapto-2-pentanone | 0.079 | 0.129 | 0.127 | 0.308 | 0.308 | 0.598 | 0.554 | 0.887 | 1.098 | 1.102 | 0.973 | 1.129 | 0.997 | 1.428 | 0.123 |
| 2-Thiophenemethanethiol | 0 | 0 | 0 | 0.004 | 0 | 0.013 | 0 | 0.098 | 0 | 0.129 | 0.119 | 0.146 | 0.114 | 0.217 | 0.079 |
| 2-Methyl-3-[(2-methyl-3-thienyl)dithio]furan | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.079 | 0.021 | 0.178 | 0.035 | 0.218 | 0.015 |
| Bis(2-furfuryl)sulfide | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.114 | 0.034 | 0.215 | 0.179 | 0.449 | 0.072 |
| **Furfural** | 0.214 | 0.428 | 0.307 | 1.106 | 0.351 | 2.976 | 2.891 | 4.028 | 3.315 | 4.894 | 3.381 | 4.128 | 3.381 | 4.683 | 5.566 |
| 1,2,3-Trithiolane | 0 | 0 | 0 | 0 | 0 | 0.336 | 0 | 0.599 | 0.156 | 0.819 | 0.208 | 1.158 | 0.179 | 1.298 | 0.200 |
| 2-Acetylthiazole | 0.477 | 0.728 | 0.634 | 0.998 | 0.708 | 1.794 | 0.814 | 3.084 | 2.488 | 3.477 | 3.079 | 4.024 | 3.398 | 4.798 | 1.299 |

Text anchors, verbatim: totals rose "from 1.023 μg/L at 20 min to 17.359 μg/L at 120 min, among which
the sulfur-containing compounds increased from 0.809 μg/L to 13.978 μg/L"; and the two thiols
"reached 1.237 μg/L and 3.734 μg/L at 120 min".

## 3. Table 1 — the three-temperature slice at 120 min, with standard deviations (µg/L)

| compound | 100 °C TTCA | 100 °C TTCA+Xyl | 120 °C TTCA | 120 °C TTCA+Xyl | 140 °C TTCA | 140 °C TTCA+Xyl |
|---|---:|---:|---:|---:|---:|---:|
| **2-Methyl-3-furanthiol** | 1.237 ± 0.142 | 1.729 ± 0.152 | 1.388 ± 0.256 | 2.498 ± 0.13 | **5.907 ± 0.085** | 3.238 ± 0.735 |
| **2-Furfurylthiol** | 3.734 ± 0.085 | 4.736 ± 0.639 | 4.107 ± 0.137 | 6.398 ± 0.266 | **11.439 ± 0.265** | 6.123 ± 0.235 |
| 3-Thiophenethiol | – | – | 0.178 ± 0.027 | 7.389 ± 0.289 | 0.268 ± 0.062 | 11.497 ± 0.649 |
| 2-Thiophenemethanethiol | 0.119 ± 0.074 | 0.146 ± 0.023 | 0.196 ± 0.015 | 0.779 ± 0.152 | 0.332 ± 0.21 | 2.479 ± 0.182 |
| 3-Mercapto-2-butanone | 1.303 ± 0.153 | 1.673 ± 0.115 | 1.798 ± 0.127 | 3.769 ± 0.274 | 2.398 ± 0.141 | 3.429 ± 0.099 |
| 3-Mercapto-2-pentanone | 0.973 ± 0.099 | 1.129 ± 0.007 | 1.298 ± 0.253 | 2.286 ± 0.142 | 1.736 ± 0.234 | 2.103 ± 0.208 |
| 2-Methyl-3-[(2-methyl-3-thienyl)dithio]furan | 0.021 ± 0.004 | 0.178 ± 0.039 | 0.098 ± 0.033 | 0.178 ± 0.031 | 0.137 ± 0.036 | 0.217 ± 0.080 |
| Bis(2-furfuryl)sulfide | 0.034 ± 0.007 | 0.215 ± 0.083 | 0.042 ± 0.004 | 0.237 ± 0.113 | 0.078 ± 0.018 | 0.339 ± 0.158 |
| **Furfural** | 3.381 ± 0.089 | 4.128 ± 0.139 | 5.793 ± 0.422 | 8.497 ± 0.514 | 11.039 ± 0.302 | 26 ± 0.725 |
| Pyrazine | – | – | – | 0.348 ± 0.061 | – | 1.389 ± 0.071 |
| Methylpyrazine | – | – | – | 1.039 ± 0.059 | – | 2.315 ± 0.129 |
| 2,5-Dimethylpyrazine | – | – | – | 2.578 ± 0.106 | – | 6.475 ± 0.741 |
| **sulfur subtotal** | 13.978 ± 0.568 | 20.19 ± 0.514 | 35.866 ± 1.382 | 91.807 ± 1.67 | 60.4 ± 1.62 | 163.432 ± 3.102 |
| **total** | 17.359 ± 0.527 | 24.318 ± 0.517 | 42.024 ± 1.737 | 106.116 ± 2.239 | 73.157 ± 1.557 | 205.954 ± 2.4 |

The pyrazine block is the sharpest structural result: **pyrazines appear only when extra xylose is
present**, at every temperature, and never in TTCA alone.

## 4. The finding that bears directly on the repository's thiol-sink question

**2-methyl-3-furanthiol rises monotonically with temperature in TTCA alone — 1.237 → 1.388 → 5.907
µg/L at 100 → 120 → 140 °C — and FALLS when xylose is added at 140 °C (3.238 against 5.907).** The
same inversion happens to 2-furfurylthiol (11.439 alone against 6.123 with xylose). Wave B27 gated
its dicarbonyl-redox couple on the premise that the disulfide deficit in fed pots is an oxidant
question; here a pot given *more* sugar makes *less* free thiol at the highest temperature, with the
difference appearing in thiophenes and mixed disulfides. That is a branching result and a candidate
premise check.

## 5. The three limits

1. **The quantification sets every response factor to 1.** Verbatim: "the detected compounds were
   quantified by the comparison of peak areas and the correction factor was 1." This is an
   internal-standard-normalised peak-area ratio in µg/L — one tier above a bare peak area, and not a
   calibrated concentration. The repository must classify it as such before scoring it.
2. **This is TTCA, not xylose + cysteine.** The only true sugar + amino-acid pot is the single MRPs
   control column. The model can charge TTCA, so this is usable — but not as a xylose–cysteine row.
3. **None of the three canonical dimers is here.** There is no bis(2-methyl-3-furyl) disulfide, no
   difurfuryl disulfide and no mixed 2-furfuryl 2-methyl-3-furyl disulfide. "Bis(2-furfuryl)sulfide"
   is a **mono**sulfide and must not be mapped onto the dimer lane. Only `wang2026b_extraction.md`
   carries those three, and only as odour-activity values.

Fig. 2 (3-deoxy- and 1-deoxyosone time courses) and Fig. 3 (regenerative cysteine, diacetyl,
methylglyoxal and glyoxal against time at all three temperatures, in mmol/L) are **figure-only** and
are the single most valuable unread object left in this paper: they are fed-intermediate dicarbonyl
time courses. Table S1 and S2 are supplementary and not on disk.
