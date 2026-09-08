# Huang 2017 — EXTRACTION (glucose or maltose + leucine or isoleucine, 10:1 molar, 0.1 mol/L phosphate pH 5.2, 90-130 C, 0-360 min; 2- and 3-methylbutanal by purge-and-trap GC-FID; a three-step multiresponse fit)
### The only paper in the corpus that fits leucine and isoleucine loss and the methylbutanals over a five-temperature ladder in water; it prints rate constants at 100 C and two barriers per system, but the absolute reactant concentrations are never stated and the printed Arrhenius lines do not reproduce the printed 100 C constants.

**Source on disk:** `data/articles/huang2016.pdf` (11 pp., owner's download, 2026-09-08; the file name says
2016 because the article was accepted February 2016 and issued in J. Food Process Eng. 2017). Read from
the text layer (`scratchpad/articles/huang2016.txt`), which came through clean for Tables 1 and 2 apart
from the exponent signs (the layer prints `5.68E 203` for 5.68E-03 and `e22953.8/T` for e^(-2953.8/T));
each of those was checked against the printed Ea (slope x R) below, so the signs are not in doubt. The
whole PDF text was searched again with `pypdf` for any absolute concentration (mmol, mol/L, g/L, mg/L):
there is none. Figures 1-7 are FIGURE-ONLY; their axis units are not in the text layer. No value was read
off any figure.

## 0. Identity

| field | value |
|---|---|
| Title | "A kinetic study on the formation of 2- and 3-methylbutanal" |
| Authors | Yarong Huang*, Johannes Tippmann, Thomas Becker (Institute of Brewing and Beverage Technology, TU München, Freising) |
| Venue | Journal of Food Process Engineering (Wiley); received 14 Aug 2015, accepted 3 Feb 2016; issued 2017 (article e12375 by the DOI; the PDF header still reads "00 (2016) 00-00") |
| DOI | 10.1111/jfpe.12375 |
| Naming | 2-MB = 2-methylbutanal (from isoleucine); 3-MB = 3-methylbutanal (from leucine); GL / ML / GI / MI = glucose-leucine / maltose-leucine / glucose-isoleucine / maltose-isoleucine; ZP = the lumped "intermediate" pool of the model (Schiff base, Amadori products, dicarbonyls); AP = "degradation products" of the aldehyde; TBA = thiobarbituric acid index (wort thermal load) |
| Text-layer quirks | `8C` = °C; `5` in equations = "="; `2` before a number = minus sign; `6` between numbers = ±; `/C1` = a multiplication dot; `lL` = µL |

## 1. Why it matters

Programme 6 (roadmap 5c) needs per-amino-acid Strecker rates or yields at two or more temperatures
in water. This paper is the corpus's only aqueous leucine and isoleucine ladder: five temperatures
(90-130 C), two sugars, six hours, with the amino acid, the sugar and the aldehyde all measured and
a rate constant per step printed at 100 C plus an Arrhenius line per step. It is the natural
temperature counterpart to Balagiannis 2010's 130 C liver-extract model and the wort-side match for
the trunk's pH 5-6 pots. What it cannot give is an absolute rate: the amino-acid concentration is
never printed (only "mol ratio 10:1"), and k1 is a pseudo-first-order loss constant that lumps every
sink of the amino acid, not the Strecker step alone.

## 2. Methods as they matter to a model

- **Reactants.** "Disaccharide sugar (maltose) or monosaccharide sugar (glucose) and amino acid
  (leucine or isoleucine), mol ratio 10:1, were dissolved in a phosphate buffer (0.1 mol/L; pH 5.2)".
  The absolute concentrations of sugar and amino acid are NOT stated anywhere in the paper (text,
  tables, captions all checked). The ratio was chosen because "the molar ratio between sugar and
  amino compound of wort is more than 3:1 (Fox et al. 1983)". Reagents Sigma-Aldrich; Milli-Q water.
- **Buffer / pH.** 0.1 mol/L phosphate, pH 5.2 (wort pH). Whether pH was followed during heating is
  not stated.
- **Vessel and heating.** "closed glass tubes" in a "heating block"; "various times (0-360 min) at
  different temperatures (90-130 C)"; tubes cooled in ice water at the planned time. From the text
  the temperatures used are 90, 100, 110, 120 and 130 C (the results discuss 90, 100, "110 and 130",
  120 and 130 C). Tube volume and headspace not stated.
- **Replicates.** "The tests were repeated three times." No error bars are printed except the ±
  on the Table 1 constants.
- **Analytics.** Sugars by HPAEC-PAD (Dionex ICS-1000, CarboPac PA10; 100 µL sample + 900 µL
  internal-standard solution; the internal standard is not named). Amino acids by OPA/FMOC
  derivatisation and HPLC (Dionex UltiMate 3000), samples diluted 1:10 by weight. 2-MB and 3-MB by
  purge-and-trap GC-FID (HP 5890, Chrompack PTI, HP-Innowax 50 m x 0.20 mm x 0.40 µm and HP Ultra-2
  50 m x 0.20 mm x 0.33 µm, purge vessel 50 C, hydrogen carrier, FID 250 C; MEBAK method after
  Pfenninger 1993 "with major modification"). **No calibration, internal standard, LOD or recovery is
  given for the aldehydes**, and the concentration unit of the aldehyde figures is not in the text.
- **Model (Scheme 3; Eqs. 9-12, verbatim structure).**
  d[Leu]/dt = -k1 [Leu];  d[ZP]/dt = k1 [Leu] - k2 [ZP];  d[3MB]/dt = k2 [ZP] - k3 [3MB];
  d[AP]/dt = k3 [3MB]. First order in each species; the sugar does not appear (it is in ten-fold
  excess and its loss is described separately as "linear with heating time"). Fitted by multiresponse
  regression in Athena Visual Studio. k3 "is the smallest and has the largest corresponding interval"
  and is NOT tabulated.
- **Statistics.** Scheffé F-test between sugar/amino-acid pairs, p < 0.05.

## 3. Tables re-typed

### Table 1. "The reaction rate constants of the reaction between GL, ML, GI and MI at 100 C"

Units as printed: "/min" (min^-1). Values are mean ± (the paper does not say whether the ± is a
standard deviation, a standard error or a 95 % interval from the regression).

| constant | GL | ML | GI | MI |
|---|---|---|---|---|
| k1 (/min) | 5.68E-03 ± 6.74E-04 | 7.88E-04 ± 1.39E-05 | 1.99E-04 ± 2.03E-05 | 3.24E-04 ± 9.82E-05 |
| k2 (/min) | 1.14E-03 ± 9.09E-05 | 1.62E-03 ± 6.76E-04 | 3.61E-04 ± 7.90E-05 | 1.20E-03 ± 2.20E-04 |

Unit reconciliation: first-order constants; 1 min^-1 = 1.667e-2 s^-1. In s^-1: k1 = 9.47e-5 (GL),
1.31e-5 (ML), 3.32e-6 (GI), 5.40e-6 (MI); k2 = 1.90e-5 (GL), 2.70e-5 (ML), 6.02e-6 (GI), 2.00e-5 (MI).
Half-lives of the amino acid at 100 C from k1: 122 min (GL), 880 min (ML), 3480 min (GI), 2140 min
(MI). The text's qualitative statements ("at 130 C more than 60 % of leucine ... degraded after 360
min", "at 90 C about 20 %") are consistent with constants of this size.

### Table 2. "Kinetic parameters with temperature dependence for the model in Schema 3"

The text layer prints the Arrhenius functions with the signs lost; they are re-typed here in the form
k = A exp(-B/T) with T in kelvin, k in min^-1 (the unit of Table 1). The "Arrhenius-function 1" row
is step 1 (k1) and "Arrhenius-function 2" is step 2 (k2); R²1 and R²2 belong to those lines.

| row | GL | ML | GI | MI |
|---|---|---|---|---|
| Arrhenius-function 1 (k1) | k = 3.015 · e^(-2953.8/T) | k = 2.036 × 10^3 · e^(-5448.8/T) | k = 3.895 × 10^4 · e^(-7000.8/T) | k = 1.622 × 10^3 · e^(-5732.6/T) |
| Arrhenius-function 2 (k2) | k = 1.434 × 10^9 · e^(-10020/T) | k = 6.392 × 10^12 · e^(-13283/T) | k = 1.892 × 10^14 · e^(-14515/T) | k = 2.651 × 10^12 · e^(-13259/T) |
| R²1 | 0.9891 | 0.9514 | 0.9355 | 0.9816 |
| R²2 | 0.9874 | 0.9664 | 0.9369 | 0.9704 |
| Ea1 (kJ/mol) | 24.56 | 45.3 | 58.2 | 47.67 |
| Ea2 (kJ/mol) | 83.31 | 110.43 | 120.68 | 110.24 |
| Etotal (kJ/mol) | 107.87 | 155.73 | 178.88 | 157.91 |

**Arithmetic checks (mine).**
- Slope × R: 2953.8 × 8.314 = 24.56; 5448.8 × 8.314 = 45.30; 7000.8 × 8.314 = 58.20; 5732.6 × 8.314
  = 47.66; 10020 × 8.314 = 83.31; 13283 × 8.314 = 110.43; 14515 × 8.314 = 120.68; 13259 × 8.314 =
  110.24 kJ/mol. Every printed Ea is the slope of its printed line, which fixes the lost signs.
- "Etotal" = Ea1 + Ea2 exactly (24.56 + 83.31 = 107.87, etc.). The sum of two consecutive first-order
  barriers is not the barrier of any measurable quantity; the abstract's "total activation energies
  ... 107.87 to 178.88 kJ/mol" should not be carried as a barrier.
- **The printed lines do not reproduce Table 1 at 100 C (373.15 K).** Evaluating each line at 373.15
  K and dividing the Table 1 value by it:

| constant | GL | ML | GI | MI |
|---|---|---|---|---|
| k1 from line at 100 C (min^-1) | 1.10e-3 | 9.27e-4 | 2.77e-4 | 3.45e-4 |
| Table 1 k1 / line | **5.16** | 0.85 | 0.72 | 0.94 |
| k2 from line at 100 C (min^-1) | 3.12e-3 | 2.22e-3 | 2.42e-3 | 9.81e-4 |
| Table 1 k2 / line | **0.36** | 0.73 | **0.15** | 1.22 |

  ML and MI agree within 30 %; GL's k1 is five times above its own Arrhenius line (the line gives
  5.68e-3 min^-1 only at 198 C) and GI's k2 is seven times below (its line gives 3.61e-4 only at 83
  C). The two tables cannot both be right for the glucose systems; Table 1 may be a single-temperature
  fit and Table 2 a global fit, but the paper does not say. Any use must carry both numbers and the
  discrepancy.
- Lines evaluated over the ladder (min^-1), for the record: GL k1 8.85e-4 (90 C) → 1.98e-3 (130 C);
  GL k2 1.49e-3 → 2.30e-2; ML k1 6.20e-4 → 2.75e-3; ML k2 8.32e-4 → 3.14e-2; GI k1 1.65e-4 → 1.12e-3;
  GI k2 8.29e-4 → 4.37e-2; MI k1 2.26e-4 → 1.08e-3; MI k2 3.69e-4 → 1.38e-2.

### Numbers in the text (the curves are FIGURE-ONLY)

- Leucine (Fig. 1): with glucose, > 60 % degraded after 360 min at 130 C, about 20 % at 90 C.
- Isoleucine (Fig. 2): with maltose, 35 % degraded after 360 min at 130 C, about 10 % at 90 C,
  "constant at the beginning or even increased slightly at 90 C"; with glucose at 110-130 C the loss
  "seem[s] to occur in three steps" (fast, slow, fast again).
- "Compared to isoleucine, leucine is degraded about two times faster than isoleucine."
- Glucose (Fig. 3): about 40 % lost after 360 min at 130 C with leucine, about 30 % with isoleucine;
  no significant difference between the two amino acids at 90-120 C; mannose, fructose, maltose,
  maltotriose and sucrose detected in the glucose samples. Maltose (Fig. 4): linear loss, "a little
  bit slower" than glucose, significantly faster with isoleucine than with leucine.
- 3-MB (Fig. 5): raising the temperature from 90 to 130 C increases the concentration "about 100
  times"; 90 and 100 C curves close, the 100 C values significantly higher. Ordering of aldehyde
  formation "GL > ML > GI > MI"; "the type of sugar ... influences the amount of Strecker flavors
  greater than the type of amino acids."
- Fig. 7 compares model and data for leucine, the intermediate (theoretical only) and 3-MB at 100
  and 120 C; no numbers.
- Secondary quotations (not this paper's data): Cremer & Eichner 2000 zero-order 10-120 min, Ea 124
  (2-MB) and 120 (3-MB) kJ/mol; Chan & Reineccius 1994 pseudo-zero order, Ea(3-MB) 80.4 kJ/mol;
  Balagiannis 2009 Ea of the first two steps 137 ± 15.2 and 48.7 ± 8.4 kJ/mol; Gomyo 1989 blue
  pigment Ea 73 kJ/mol. Pointers only.

## 4. Kinetic numbers the repository can use

Registry (`data/keys/compounds.yml`): 3-methylbutanal → `3_methylbutanal`; 2-methylbutanal →
`2_methylbutanal`; leucine, isoleucine, glucose, maltose, the lumped intermediate → not in registry
(reaction_rules.yml uses the short names Leu, Ile, Glc).

| quantity | value | unit | conditions | reaction order (authors) | source location | evidence class |
|---|---|---|---|---|---|---|
| k1, leucine loss with glucose (GL) | 5.68e-3 ± 6.74e-4 | min^-1 | 0.1 mol/L phosphate pH 5.2, sugar:amino acid 10:1 (absolute concentrations not printed), 100 C, closed tubes, 0-360 min, n = 3 | pseudo-first order in the amino acid | Table 1 | measured_rate (lumped amino-acid loss, all sinks) |
| k1, ML / GI / MI | 7.88e-4 ± 1.39e-5 / 1.99e-4 ± 2.03e-5 / 3.24e-4 ± 9.82e-5 | min^-1 | same | same | Table 1 | measured_rate (lumped) |
| k2, intermediate → 3-MB (GL) / → 3-MB (ML) / → 2-MB (GI) / → 2-MB (MI) | 1.14e-3 ± 9.09e-5 / 1.62e-3 ± 6.76e-4 / 3.61e-4 ± 7.90e-5 / 1.20e-3 ± 2.20e-4 | min^-1 | same, 100 C | first order in an unmeasured pool | Table 1 | measured_rate (model-identified only; see Flag 3) |
| Ea1 (k1), GL / ML / GI / MI | 24.56 / 45.3 / 58.2 / 47.67 | kJ/mol | 90-130 C, five temperatures | Arrhenius on the fitted k1 | Table 2 | measured_barrier (lumped loss; R² 0.94-0.99) |
| Ea2 (k2), GL / ML / GI / MI | 83.31 / 110.43 / 120.68 / 110.24 | kJ/mol | same | Arrhenius on the fitted k2 | Table 2 | measured_barrier (see Flag 2: the lines and Table 1 disagree for GL and GI) |
| Arrhenius prefactors A1, A2 | 3.015, 2.036e3, 3.895e4, 1.622e3 (k1); 1.434e9, 6.392e12, 1.892e14, 2.651e12 (k2) | min^-1 | same | — | Table 2 | measured (printed with the lines) |
| k1(GL)/k1(GI) at 100 C; k1(ML)/k1(MI) | 28.5; 2.43 | — | same | — | derived from Table 1 | within_study_ratio (leucine vs isoleucine; the text says "about two times", see Flag 4) |
| k1(GL)/k1(ML); k1(GI)/k1(MI) | 7.2; 0.61 | — | same | — | derived from Table 1 | within_study_ratio (glucose vs maltose; the sign reverses between amino acids) |
| aldehyde formation ordering | GL > ML > GI > MI | — | all temperatures | — | text (Results, "2-MB and 3-MB") | within_study_ratio (ordinal only) |
| 3-MB level at 130 C vs 90 C, same time | "about 100 times" | — | GL | — | text | within_study_ratio (approximate, figure-based statement by the authors) |
| leucine loss at 360 min, 130 C / 90 C (GL) | > 60 % / about 20 % | % | GL | — | text | level_only (approximate) |
| isoleucine loss at 360 min, 130 C / 90 C (MI) | 35 % / about 10 % | % | MI | — | text | level_only |
| glucose loss at 360 min, 130 C | about 40 % (GL), about 30 % (GI) | % | — | — | text | level_only |
| k3 (aldehyde loss) | not tabulated ("the smallest ... largest interval") | — | — | first order | text | not reported |
| 2-MB, 3-MB, amino acid, sugar concentrations vs time | — | unit not in text | 90-130 C | — | Figs. 1-7 | figure_only |

Cross-reference inside the repo: Balagiannis 2010 (`balagiannis2010_extraction.md`) fits the same
Leu/Ile → intermediate → methylbutanal chain in a liver extract at 130 C, and Parker 2013 quotes the
2009 companion (JAFC 57:9916, not on disk) at 120/130/140 C. `jousse2002_extraction.md` carries the
lumped Strecker step. The rule R07 (`data/lit/reaction_rules.yml`) is the mechanistic step; this
paper's k1 is upstream of it (amino acid loss into the Amadori/dicarbonyl pool), and k2 is the
lumped R07 + release step.

## 5. Flags

1. **No absolute concentrations.** Only "mol ratio 10:1" and "0.1 mol/L phosphate" are printed. The
   pseudo-first-order k1 is usable as a rate per unit amino acid only if the ten-fold sugar excess is
   accepted as constant; it cannot be turned into a second-order constant (dicarbonyl × amino acid,
   the B18 convention) without the sugar concentration, and the aldehyde curves cannot be converted to
   mmol/L. Yields relative to the amino acid are recoverable from the figures only, which the
   repository does not transcribe.
2. **Table 1 and Table 2 disagree for the glucose systems.** The printed Arrhenius line for GL gives
   k1 = 1.10e-3 min^-1 at 100 C against Table 1's 5.68e-3 (factor 5.2), and GI's k2 line gives
   2.42e-3 against 3.61e-4 (factor 6.7); ML and MI agree within 30 %. The paper does not say how the
   two tables were obtained (per-temperature vs global fit). Neither table should be used alone.
3. **k1 is not a Strecker rate.** It is the total loss of the amino acid: Amadori formation, Strecker,
   melanoidin binding. With a 10:1 sugar excess at pH 5.2 most of the amino acid consumed is not
   released as aldehyde (the aldehyde yields of Figs. 5-6 are not printed, but the authors' own
   Fig. 7 shows the "intermediate" pool accumulating). k2 is identified only through the model
   (the intermediate was never measured), so its value and its Ea depend on the assumed chain. The
   low Ea1 of GL (24.6 kJ/mol) is what a diffusion- or equilibrium-limited lumped step looks like,
   not an activation barrier of a bond-breaking step.
4. **Leucine/isoleucine ratio: 28.5 (Table 1, glucose) vs "about two times" (text) vs 2.4 (maltose).**
   The text's factor two matches the maltose pair only. Any within-study ratio taken from this paper
   should be the maltose pair or should carry both.
5. **Aldehyde quantification undocumented.** No internal standard, calibration, LOD or recovery for
   2-MB/3-MB; purge-and-trap at 50 C with FID; the concentration unit of Figs. 5-6 is not in the
   text. The "about 100 times" between 90 and 130 C is the authors' reading of their figure.
6. **The model has no sugar term and no dicarbonyl.** Rate constants are for the pH 5.2, 10:1
   regime only; the Martins-type dependence on initial concentrations (which the introduction
   itself cites as "no effect") is untested here.
7. **Open tubes' volatility.** 2-MB/3-MB boil at ~91 C; the tubes were closed, but headspace volume
   is not stated, so the liquid-phase aldehyde at 100-130 C depends on the unknown headspace ratio.
8. **k3 not reported**, so the aldehyde-loss step (which Balagiannis 2010 needed) has no number here.
9. **"Etotal"** is an arithmetic sum of two step barriers, not a measurable barrier; the abstract's
   107.87-178.88 kJ/mol range must not be carried as an Ea of methylbutanal formation.
10. Registry gap: leucine, isoleucine and the sugars have no `compounds.yml` rows (the rules file
    uses short names); the two aldehydes are registered.
