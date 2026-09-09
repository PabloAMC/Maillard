# Leahy & Reineccius 1989 (ch. 7) — EXTRACTION (0.1 M sugar + 0.1 M amino acid, pH 9.0 borate, 75/85/95 C, up to 24 h; lysine or asparagine x glucose, fructose or ribose)
### Part One of the pair: pseudo-zero-order pyrazine formation rates in ppm/h and 3-point Arrhenius Ea for five sugar-amino acid systems.

**Source on disk:** `data/articles/leahy1989.pdf` (16 pp., owner's download, 2026-09-08). ACS Symposium
Series 388 (*Flavor Chemistry: Trends and Developments*, Teranishi, Buttery, Shahidi eds.), chapter 7,
pp. 76-91. The text layer is an OCR layer with glyph errors in prose ("diirethylpyrazine",
"œncentration") but the numbers in Tables I and II came through cleanly; Table III's columns were
scrambled in the text layer. **Tables I, II and III were verified cell by cell against 130-dpi rasters
of printed pages 83, 84 and 87** (`scratchpad/img/l89-08.png`, `l89-09.png`, `l89-12.png`); every
number below matches the raster. Figures 2 and 3 (concentration-time plot, distribution bars) are
FIGURE-ONLY and were not read.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics of Formation of Alkylpyrazines. Effect of Type of Amino Acid and Type of Sugar" |
| Authors | M. M. Leahy (Ocean Spray Cranberries, Middleboro MA) and G. A. Reineccius (Univ. of Minnesota, St. Paul) |
| Venue | ACS Symposium Series 388, *Flavor Chemistry: Trends and Developments*, 1989, ch. 7, pp. 76-91. Received September 23 1988; published February 21 1989 |
| DOI | 10.1021/bk-1989-0388.ch007 |
| Naming | "2-methylpyrazine" = methylpyrazine; "lys-glu", "lys-fruc", "lys-rib", "asp-glu", "asp-fruc" in Table II mean lysine/asparagine with glucose/fructose/ribose ("asp" is asparagine, not aspartic acid); Table III uses L/A for lysine/asparagine and G/F/R for the sugars |
| Companion | ch. 18 of ACS Symp. Ser. 409 (`leahy1989a_extraction.md`): same lysine-glucose pH 9.0 data re-printed as its pH 9.0 arm, plus pH 5.0 / 7.0 and a water-activity study. Balagiannis 2015 cites this chapter as "Leahy & Reineccius 1989a" |
| Cited in repo | `apriyantono1993_extraction.md` and `k5b_dmhf_synthesis.md` name the ch. 18 chapter as a fetch target; `chan1994_extraction.md` records the pyrazine Ea-vs-pH claim second-hand and says "go to the primaries" |

## 1. Why it matters

This is one of very few sources with **printed pyrazine formation rates (slope, intercept, n, r2) at
three temperatures and an Arrhenius Ea for each compound**, for five sugar-amino acid pairs, all
under one analytical method with an internal standard and empirically determined response factors.
For a pyrazine lane it supplies: (i) the temperature dependence (Ea 114-187 kJ/mol, i.e. 27-45
kcal/mol) of pseudo-zero-order pyrazine, methylpyrazine and dimethylpyrazine formation from
sugar + amino acid; (ii) within-study sugar and amino-acid ratios; (iii) the 2 h / 95 C product
distribution. Its limits are equally clear: pH 9.0 borate, 75-95 C (below the engine's 100-145 C
window), lysine and asparagine only (no glycine, no alanine), conversions of order 0.1 %, and only
three temperatures per Ea. Nothing here isolates the Strecker or condensation step: the rates are
whole-cascade rates from sugar + amino acid.

## 2. Methods as they matter to a model

- **Reactants.** "asparagine-fructose, asparagine-glucose, lysine-fructose, lysine-glucose and
  lysine-ribose at concentrations of 0.1M for both amino acid and sugar, in a pH 9.0 0.1M borate
  buffer". So **100 mmol/L sugar + 100 mmol/L amino acid, 1:1**, 0.1 M borate. (Ch. 18 specifies
  L-lysine monohydrochloride for the lysine runs; this chapter says only "lysine".) A cysteine +
  glucose control (0.1 M each, pH 9.0, 95 C, 2 h) gave no detectable pyrazines.
- **Vessel / heating.** 10 mL in Teflon-capped 25 mm o.d. x 150 mm Pyrex test tubes, water bath at
  75, 85 and 95 C, "for up to 24 h". No stirring, air headspace, not sealed against volatile loss
  beyond the cap.
- **Sampling.** "Samples were taken at 7 to 8 time intervals. Eighteen to 20 total samples per
  temperature were analyzed." "Two zero points were used for each regression. Duplicate samples were
  tested at the early sampling times vs. triplicate samples at later times, as variations in
  concentration among replicates increased with increased reaction time. Each data point collected
  was treated separately in the regression analyses." The actual sampling times are not printed
  (Figure 2 has them on its axis; figure-only).
- **pH.** "Although the amino acid/sugar solutions were buffered ... a drop in pH was encountered
  with increasing reaction times" — magnitude not printed. After heating each sample was
  re-adjusted to pH 9.0 with 0.1 N NaOH before analysis (an analytical step; the reaction pH
  drifted downward from 9.0 by an unstated amount).
- **Internal standard.** 1 mL of 2 ppm 2-methoxypyrazine in water (= 2 µg) added after heating;
  final sample volume 15 mL (10 mL reaction + 1 mL IS + NaOH + water to 15 mL).
- **Isolation / GC.** Automated purge and trap (HP 7675A): 15 mL sample purged 10 min with H2 at 90
  mL/min onto a 4" x 1/4" Tenax precolumn; desorbed 3 min at 180 C, split 50:1, cryofocused in
  liquid N2; DB-225 25 m x 0.32 mm i.d.; isothermal 50 C then 200 C for 1.5 min; HP 5880A with
  **nitrogen-phosphorus detector**; injector 225 C, NPD 280 C. Run time about 30 min.
- **Quantification: absolute, internal-standard method with empirical response factors.**
  Amt_C = Amt_ISTD x AC_C / (RF x AC_ISTD), amounts in µg/mL, RF "empirically determined by adding
  known amounts of each compound to 15 ml of pH 9.0 borate buffer, followed by a purge of the sample
  under usual conditions of analysis" — i.e. the calibration includes the purge-and-trap recovery in
  the same matrix minus reactants. Reported unit **ppm = µg/mL**. Whether "ppm" refers to the 15 mL
  analysed volume or has been corrected back to the 10 mL reaction volume is not stated (factor
  1.5; see Flags). LOD not stated; the text reports "3 ppb" of 2,3-dimethylpyrazine in the
  asparagine-glucose system, so the working LOD is at or below 3 µg/L.
- **Identification.** Co-chromatography with authentic standards (Pyrazine Specialties), confirmed by
  GC-MS (Carlo Erba / Kratos MS 25, 70 eV) on a 0.1 M glucose-lysine sample heated 6 h at 95 C.
  Compounds resolved (Figure 1): pyrazine, 2-methylpyrazine, 2,5-dimethylpyrazine,
  2,6-dimethylpyrazine, 2,3-dimethylpyrazine (+ the IS). No ethyl- or trimethylpyrazines reported.
- **Kinetic treatment (verbatim).** "dA/dt = k A^n ... A = A0 + kt for a zero order reaction ...
  The formation of pyrazines appears to better fit a pseudo zero order reaction rather than first
  order reaction. Plotting concentrations of pyrazines formed versus time of reaction gave the better
  fit of the line, usually with a coefficient of determination (r2) of greater than 0.95. For a
  pseudo first order reaction, a curve rather than a line was obtained. General least squares
  analysis of the data was used to compute rate constants." Units: A in ppm, t in h, so **k in
  ppm/h**. Arrhenius: k = k0 exp(-Ea/RT) with R = 1.986 cal/mol/K, **Ea in kcal/mol**, three
  temperatures, so three points per Ea. The authors' own caveat: "since three points were used in
  the Arrhenius plots ... further research in this area is necessary to determine the significance
  of these findings", and the ANOVA on Ea "assumption [of constant variance] was found to be
  violated".

## 3. Tables re-typed

### Table I. "Regressions for the effect of type of sugar and amino acid on the formation of pyrazines"

Columns as printed: model system, temperature, k (ppm/hr), k0 (intercept), number of samples, r2.
The last column is added here: **k converted to µmol L-1 min-1 = k(ppm/h) / MW x 1000 / 60**, with
MW pyrazine 80.088, methylpyrazine 94.115, dimethylpyrazines 108.141 g/mol (1 ppm = 1 mg/L). The
conversion assumes "ppm" is the reaction-mixture concentration (Flag 3).

| system | compound | T (C) | k (ppm/h) | intercept (ppm) | n | r2 | k (µmol L-1 min-1) |
|---|---|---:|---:|---:|---:|---:|---:|
| lysine-glucose | pyrazine | 95 | 3.596 | 0.0596 | 22 | 0.994 | 0.748 |
| | | 85 | 0.490 | 0.458 | 22 | 0.960 | 0.102 |
| | | 75 | 0.214 | 0.279 | 20 | 0.965 | 0.0445 |
| | 2-methylpyrazine | 95 | 2.837 | -0.104 | 22 | 0.995 | 0.502 |
| | | 85 | 0.422 | 0.142 | 22 | 0.967 | 0.0747 |
| | | 75 | 0.159 | 0.0910 | 20 | 0.941 | 0.0282 |
| | 2,5-dimethylpyrazine | 95 | 0.186 | -0.0457 | 20 | 0.995 | 0.0287 |
| | | 85 | 0.0247 | -0.00604 | 22 | 0.985 | 0.00381 |
| | | 75 | 0.00668 | -0.00536 | 16 | 0.942 | 0.00103 |
| | 2,3-dimethylpyrazine | 95 | 0.0229 | -0.00569 | 16 | 0.948 | 0.00353 |
| | | 85 | 0.00309 | -0.00173 | 18 | 0.978 | 0.000476 |
| | | 75 | 0.000677 | -0.00860 | 12 | 0.958 | 0.000104 |
| lysine-fructose | pyrazine | 95 | 1.359 | -0.226 | 20 | 0.995 | 0.283 |
| | | 85 | 0.395 | -0.0991 | 22 | 0.995 | 0.0822 |
| | | 75 | 0.134 | -0.103 | 21 | 0.963 | 0.0279 |
| | 2-methylpyrazine | 95 | 1.105 | 0.301 | 20 | 0.995 | 0.196 |
| | | 85 | 0.388 | 0.119 | 22 | 0.945 | 0.0687 |
| | | 75 | 0.116 | -0.0343 | 21 | 0.968 | 0.0205 |
| | 2,5-dimethylpyrazine | 95 | 0.175 | 0.0488 | 20 | 0.945 | 0.0270 |
| | | 85 | 0.0376 | 0.0211 | 20 | 0.951 | 0.00580 |
| | | 75 | 0.00779 | 0.00570 | 15 | 0.935 | 0.00120 |
| | 2,3-dimethylpyrazine | 95 | 0.0441 | 0.0127 | 16 | 0.954 | 0.00680 |
| | | 85 | **0.117** (sic; see Flag 1) | 0.00328 | 16 | 0.969 | (0.0180 as printed) |
| lysine-ribose | pyrazine | 95 | 3.488 | 0.530 | 21 | 0.974 | 0.726 |
| | | 85 | 0.992 | 0.764 | 22 | 0.916 | 0.206 |
| | | 75 | 0.310 | 0.747 | 22 | 0.882 | 0.0645 |
| | 2-methylpyrazine | 95 | 5.364 | 1.380 | 19 | 0.957 | 0.950 |
| | | 85 | 1.768 | 0.914 | 22 | 0.929 | 0.313 |
| | | 75 | 0.440 | 0.350 | 22 | 0.952 | 0.0779 |
| | 2,5-dimethylpyrazine | 95 | 0.152 | 0.0265 | 18 | 0.975 | 0.0234 |
| | | 85 | 0.0440 | 0.0204 | 17 | 0.933 | 0.00678 |
| | | 75 | 0.0128 | 0.00717 | 18 | 0.984 | 0.00197 |
| | 2,3-dimethylpyrazine | 95 | 0.0166 | 0.00127 | 17 | 0.973 | 0.00256 |
| | | 85 | 0.00427 | 0.000369 | 16 | 0.958 | 0.000658 |
| asparagine-glucose | pyrazine | 95 | 0.103 | -0.0371 | 22 | 0.987 | 0.0214 |
| | | 85 | 0.0422 | -0.0385 | 22 | 0.941 | 0.00878 |
| | | 75 | 0.00981 | -0.0187 | 22 | 0.947 | 0.00204 |
| | 2-methylpyrazine | 95 | 0.442 | -0.215 | 22 | 0.992 | 0.0783 |
| | | 85 | 0.179 | -0.242 | 22 | 0.926 | 0.0317 |
| | | 75 | 0.0343 | -0.102 | 22 | 0.942 | 0.00607 |
| | 2,5-dimethylpyrazine | 95 | 0.0871 | -0.0117 | 22 | 0.982 | 0.0134 |
| | | 85 | 0.0202 | -0.00255 | 22 | 0.968 | 0.00311 |
| | | 75 | 0.00455 | -0.0153 | 14 | 0.863 | 0.000701 |
| | 2,6-dimethylpyrazine | 95 | 0.0997 | -0.0875 | 21 | 0.979 | 0.0154 |
| | | 85 | 0.0188 | -0.0294 | 11 | 0.821 | 0.00290 |
| | 2,3-dimethylpyrazine | 95 | 0.00231 | -0.000192 | 20 | 0.992 | 0.000356 |
| asparagine-fructose | pyrazine | 95 | 0.0227 | -0.0135 | 20 | 0.904 | 0.00472 |
| | | 85 | 0.00653 | -0.00631 | 20 | 0.881 | 0.00136 |
| | | 75 | 0.00266 | -0.00328 | 20 | 0.947 | 0.000554 |
| | 2-methylpyrazine | 95 | 0.632 | -0.331 | 22 | 0.927 | 0.112 |
| | | 85 | 0.163 | -0.197 | 22 | 0.926 | 0.0289 |
| | | 75 | 0.0374 | -0.0928 | 22 | 0.896 | 0.00662 |
| | 2,5-dimethylpyrazine | 95 | 0.487 | -0.107 | 22 | 0.983 | 0.0751 |
| | | 85 | 0.122 | -0.0192 | 22 | 0.995 | 0.0188 |
| | | 75 | 0.0261 | -0.0323 | 20 | 0.940 | 0.00402 |
| | 2,6-dimethylpyrazine | 95 | 0.541 | -0.248 | 20 | 0.960 | 0.0834 |
| | | 85 | 0.152 | -0.104 | 20 | 0.980 | 0.0234 |
| | | 75 | 0.0356 | -0.0547 | 18 | 0.921 | 0.00549 |

Rows absent from the table (no regression printed): lysine-glucose and lysine-fructose
2,6-dimethylpyrazine (not detected in lysine systems except with ribose, per text); lysine-ribose
2,6-dimethylpyrazine and 75 C 2,3-dimethylpyrazine; asparagine-glucose 75 C 2,6-DMP and 85/75 C
2,3-DMP; asparagine-fructose 2,3-DMP (not detected).

### Table II. "Activation energies for formation of pyrazines in 0.1M sugar-amino acid systems" (Ea in kcal/mol)

kJ/mol column added (x 4.184). Last column: my unweighted 3-point Arrhenius refit of the Table I k
values, to check that Table II derives from Table I (it does, to within 0.2 kcal/mol everywhere).

| compound | system | Ea printed (kcal/mol) | Ea (kJ/mol) | refit from Table I (kJ/mol) |
|---|---|---:|---:|---:|
| pyrazine | lys-glu | 35.8 | 149.8 | 149.7 |
| | lys-fruc | 29.5 | 123.4 | 123.3 |
| | lys-rib | 30.8 | 128.9 | 128.9 |
| | asp-glu | 30.0 | 125.5 | 125.5 |
| | asp-fruc | 27.3 | 114.2 | 114.0 |
| 2-methylpyrazine | lys-glu | 36.6 | 153.1 | 153.0 |
| | lys-fruc | 28.7 | 120.1 | 120.1 |
| | lys-rib | 31.9 | 133.5 | 133.3 |
| | asp-glu | 32.6 | 136.4 | 136.5 |
| | asp-fruc | 36.0 | 150.6 | 150.7 |
| 2,5-dimethylpyrazine | lys-glu | 42.3 | 177.0 | 176.8 |
| | lys-fruc | 39.6 | 165.7 | 165.8 |
| | lys-rib | 31.5 | 131.8 | 131.8 |
| | asp-glu | 37.6 | 157.3 | 157.3 |
| | asp-fruc | 37.3 | 156.1 | 156.0 |
| 2,6-dimethylpyrazine | asp-fruc | 34.7 | 145.2 | 145.0 |
| 2,3-dimethylpyrazine | lys-glu | 44.8 | 187.4 | 187.3 |

Text: "Activation energies for alkylpyrazine formation were calculated from the slope of Arrhenius
plots, ranging from 27 to 45 kcal/mole" (table range 27.3-44.8). Duncan's test: "activation energies
for 2,5-dimethylpyrazine were significantly higher than those of pyrazine and 2-methylpyrazine",
with the variance-assumption caveat quoted in §2. Literature comparators quoted by the authors:
browning Ea 15.5 kcal/mol (glycine-glucose, Stamp & Labuza 1983) to 33 kcal/mol (IM model food,
Warmbier 1976); lysine and glucose loss 25 kcal/mol (Warmbier 1976).

**Arrhenius curvature check (mine, lysine-glucose pyrazine):** pairwise Ea 95/85 C = 218 kJ/mol,
85/75 C = 86 kJ/mol, 95/75 C = 150 kJ/mol. The 85 C point sits well below the line through the other
two; the same pattern holds for lysine-glucose methylpyrazine (209 / 101 kJ/mol). With three points
the printed Ea are averages over strongly non-linear Arrhenius plots. Use them as 75-95 C apparent
barriers only.

### Table III. "Effect of type of amino acid and type of sugar on pyrazine distributions, 2 h treatment at 95 C"

Percent of total pyrazines; last row total in ppm. Verified against the raster (the text layer
scrambled the columns). L = lysine, A = asparagine, G = glucose, F = fructose, R = ribose.

| compound | L-G | L-F | L-R | A-G | A-F |
|---|---:|---:|---:|---:|---:|
| pyrazine | 55.8 | 42.3 | 36.8 | 18.7 | 1.0 |
| 2-methylpyrazine | 41.5 | 48.6 | 61.1 | 58.8 | 29.1 |
| 2,5-dimethylpyrazine | 2.4 | 7.3 | 1.6 | 17.3 | 37.6 |
| 2,6-dimethylpyrazine | — | — | 0.3 | 4.8 | 32.3 |
| 2,3-dimethylpyrazine | 0.3 | 1.8 | 0.2 | 0.4 | — |
| **TOTAL (ppm)** | **13.1** | **5.7** | **19.9** | **0.74** | **2.2** |

Column sums: 100.0 / 100.0 / 100.0 / 100.0 / 100.0. Text: "no pyrazines were detected" for 0.1 M
cysteine + glucose, pH 9.0, 95 C, 2 h; 2,3-dimethylpyrazine "only 3 ppb in the asparagine-glucose
system versus 30-100 ppb in the lysine systems".

**Molar reconciliation of the L-G total (mine):** 13.1 ppm x 0.558 / 80.088 = 0.0913 mmol/L pyrazine;
x 0.415 / 94.115 = 0.0578 mmol/L methylpyrazine; x 0.027 / 108.141 = 0.0033 mmol/L dimethylpyrazines;
total **0.152 mmol/L pyrazines after 2 h at 95 C from 100 mmol/L lysine + 100 mmol/L glucose**, i.e.
0.30 mmol/L amino-N (two N per ring) = **0.3 % of the lysine**. Conversion is negligible over the
window, which is why zero order fits: the rates are initial rates of the whole cascade.

## 4. Kinetic numbers the repository can use

Registry mapping (`data/keys/compounds.yml`): 2-methylpyrazine -> `methylpyrazine`;
2,5-dimethylpyrazine -> `2_5_dimethylpyrazine`; 2,6-dimethylpyrazine -> `2_6_dimethylpyrazine`;
2,3-dimethylpyrazine -> `2_3_dimethylpyrazine`; **pyrazine (parent) -> not in registry as a molecule**
(class row `pyrazines` only); 2-methoxypyrazine (IS) -> class `methoxypyrazines` only; lysine,
asparagine, cysteine, glucose, fructose, ribose -> not in `compounds.yml` (`reactive_lysine` is a
marker, not the amino acid).

All rate rows share: 100 mmol/L sugar + 100 mmol/L amino acid, 0.1 M borate, initial pH 9.0
(drifting down, unquantified), capped Pyrex tube in a water bath, absolute IS quantification by
purge-and-trap NPD, k = slope of concentration vs time (pseudo-zero-order in product), 18-22 points
over up to 24 h.

| quantity | value | unit | conditions | reaction order | source location | evidence class |
|---|---|---|---|---|---|---|
| pyrazine formation rate, lysine-glucose | 3.596 / 0.490 / 0.214 (= 0.748 / 0.102 / 0.0445 µmol L-1 min-1) | ppm/h | 95 / 85 / 75 C, pH 9 | pseudo-zero (product) | Table I, p. 83 | measured_rate |
| methylpyrazine rate, lysine-glucose | 2.837 / 0.422 / 0.159 (0.502 / 0.0747 / 0.0282) | ppm/h | same | pseudo-zero | Table I | measured_rate |
| 2,5-dimethylpyrazine rate, lysine-glucose | 0.186 / 0.0247 / 0.00668 (0.0287 / 0.00381 / 0.00103) | ppm/h | same | pseudo-zero | Table I | measured_rate |
| 2,3-dimethylpyrazine rate, lysine-glucose | 0.0229 / 0.00309 / 0.000677 | ppm/h | same | pseudo-zero | Table I | measured_rate |
| all other Table I rows (lysine-fructose, lysine-ribose, asparagine-glucose, asparagine-fructose) | as re-typed in §3 | ppm/h | 95 / 85 / 75 C, pH 9 | pseudo-zero | Table I, pp. 83-84 | measured_rate (except the 0.117 cell: suspect) |
| Ea, pyrazine, lysine-glucose | 35.8 (149.8) | kcal/mol (kJ/mol) | 75-95 C, 3 points | Arrhenius on zero-order k | Table II, p. 85 | measured_barrier (3-point, curved plot) |
| Ea, methylpyrazine, lysine-glucose | 36.6 (153.1) | kcal/mol (kJ/mol) | same | same | Table II | measured_barrier |
| Ea, 2,5-dimethylpyrazine, lysine-glucose | 42.3 (177.0) | kcal/mol (kJ/mol) | same | same | Table II | measured_barrier |
| Ea, 2,3-dimethylpyrazine, lysine-glucose | 44.8 (187.4) | kcal/mol (kJ/mol) | same | same | Table II | measured_barrier |
| Ea, other 13 compound x system pairs | 27.3-39.6 kcal/mol as re-typed in §3 | kcal/mol | same | same | Table II | measured_barrier |
| sugar ratio at 95 C, pyrazine rate: ribose / glucose / fructose (lysine) | 3.488 / 3.596 / 1.359 = 0.97 : 1 : 0.38 | — | pH 9, 95 C | — | Table I | within-study ratio |
| sugar ratio at 95 C, methylpyrazine rate (lysine): ribose / glucose / fructose | 5.364 / 2.837 / 1.105 = 1.89 : 1 : 0.39 | — | same | — | Table I | within-study ratio |
| amino-acid ratio at 95 C, glucose: lysine / asparagine — pyrazine, methylpyrazine, 2,5-DMP | 3.596/0.103 = 34.9; 2.837/0.442 = 6.4; 0.186/0.0871 = 2.1 | — | same | — | Table I | within-study ratio |
| total pyrazines after 2 h at 95 C: L-G / L-F / L-R / A-G / A-F | 13.1 / 5.7 / 19.9 / 0.74 / 2.2 (L-G = 0.152 mmol/L) | ppm | pH 9 | — | Table III, p. 87 | level_only (end-of-cook; validation) |
| distribution after 2 h at 95 C | Table III percentages | % of total | pH 9 | — | Table III | within-study ratio |
| cysteine + glucose, 0.1 M each, pH 9, 95 C, 2 h | no pyrazines detected | — | — | — | text p. 85 | level_only (null) |
| 2,3-dimethylpyrazine, asparagine-glucose vs lysine systems, 2 h 95 C | 3 ppb vs 30-100 ppb | µg/L | pH 9 | — | text p. 87 | level_only |
| concentration-time points behind Table I; Figure 2 (lysine-glucose pyrazine vs time at 3 T); Figure 3 (distribution bars) | — | — | — | — | Figures 2-3 | figure_only |

Cross-reference: `data/lit/arrhenius_params.yml` `pyrazine_condensation` Ea 138.072 kJ/mol
(`estimated`, = 33.0 kcal/mol) falls inside this paper's 114-187 kJ/mol whole-cascade range; the two
are not the same step. Zhou 2024 (`zhou2024_extraction.md`) gives 100.6 / 111.7 kJ/mol for the
fed-dicarbonyl sub-cascade at pH 8, 100-120 C.

## 5. Flags

1. **Suspect cell: lysine-fructose 2,3-dimethylpyrazine at 85 C, k = 0.117 ppm/h** (raster
   confirmed as printed). It is 2.7x the 95 C value (0.0441) and would imply a negative Ea; every
   other series rises 3-6x per 10 C. Almost certainly a misprint for 0.0117 (which would give a
   2-point Ea of 146 kJ/mol, in line with its neighbours). No Ea is printed for this pair. Do not use
   the cell.
2. **The ch. 18 re-print of the lysine-glucose pH 9.0 block differs in five cells**: 2,5-DMP 95 C r2
   0.995 here vs 0.976 there; 2,3-DMP intercepts -0.00569 / -0.00173 / -0.00860 here vs -0.0057 /
   -0.00017 / -0.00086 there; Ea 2-methylpyrazine 36.6 here vs 36.7 there; Ea 2,5-DMP 42.3 here vs
   41.9 there. Rate constants k and n agree exactly. The intercept differences look like dropped
   digits in one of the two typescripts; the Ea differences (0.1 and 0.4 kcal/mol) are unexplained
   (my refit of the identical k values gives 36.6 and 42.3, i.e. this chapter's values).
3. **"ppm" volume basis unstated.** Amounts are computed in µg/mL of the 15 mL analysed sample (10 mL
   reaction mixture + IS + NaOH + water). If not corrected back, reaction-mixture concentrations and
   rates are 1.5x the printed values. The within-study ratios and all Ea are unaffected.
4. **Pseudo-zero order is an initial-rate statement.** Conversion after 2 h at 95 C is ~0.3 % of
   lysine (§3); reactants are constant, so no order in reactants is determined. Store k as a rate at
   [sugar] = [amino acid] = 0.1 M.
5. **Three-point Arrhenius with visible curvature** (pairwise 218 vs 86 kJ/mol for lysine-glucose
   pyrazine). The authors say so themselves. The Ea are 75-95 C apparent values; extrapolation to
   100-145 C is not licensed by the data.
6. **Non-zero intercepts**, some large relative to the 2 h yield (lysine-ribose methylpyrazine
   intercept 1.38 ppm at 95 C; lysine-glucose pyrazine 0.458 ppm at 85 C). Two zero points were
   forced into each regression, so a positive intercept means the early points lie above the line —
   the time course is not linear from t = 0 (fast initial phase or curvature), and the "k" is a
   secant slope over the sampled window.
7. **pH drift unquantified**; borate at pH 9.0 (borate complexes sugar diols — a matrix effect on
   sugar reactivity relative to phosphate systems; my note, not the authors'). Lysine salt form not
   stated in this chapter (monohydrochloride in ch. 18).
8. **Sampling times not printed**; "up to 24 h" and 7-8 intervals; replicates 2-3 per point but no SD
   or CI on any k or Ea; "n" ranges 11-22 because some points fell below detection (the 11-point
   asparagine-glucose 2,6-DMP 85 C row, r2 0.821, is the weakest).
9. **Prose vs Table III swap**: the text says "a greater relative yield for asparagine-fructose (49%)
   versus lysine-fructose (29%)" for 2-methylpyrazine; Table III has L-F 48.6 % and A-F 29.1 %. The
   table is self-consistent (columns sum to 100); trust the table.
10. **Text layer glyphs**: "asp" = asparagine; the OCR gives "LYSINE-GIUQOSE", "LESINE-FRUCTOSE",
    "IVSINE-RIBOSE"; all rows re-typed from the raster.
11. Not measured: sugar or amino-acid loss, any dicarbonyl, Strecker aldehydes, browning; no
    glycine or alanine system; no temperature at or above 100 C; no headspace/loss control for the
    capped tubes over 24 h at 95 C.
