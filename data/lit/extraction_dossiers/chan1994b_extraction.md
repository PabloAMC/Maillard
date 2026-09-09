# Chan 1994b — EXTRACTION (glucose 1.25 mol/L + methionine, phenylalanine, proline and leucine 0.1875 mol/L each in 400 mL of nominally 0.1 mol/L phosphate, pH 6 / 7 / 8, 75-115 C, 5 min to 7.5 h, sealed stirred Parr reactor; pseudo-zero-order formation of isovaleraldehyde, phenylacetaldehyde, 2-acetyl-1-pyrroline, 2-acetylfuran, DDMP and 5-methyl-2-phenyl-2-hexenal by GC-AED; one table of averaged activation energies, no rate constants)
### The aqueous 75-115 C Strecker-aldehyde barrier everyone quotes (isovaleraldehyde 19.2 kcal/mol = 80.3 kJ/mol, phenylacetaldehyde 21.5 kcal/mol = 90.0 kJ/mol) comes from this nine-page chapter, which prints six averaged barriers and not a single rate constant, concentration or standard deviation.

**Source on disk:** `data/articles/chan2005.pdf` (9 pp., 1,039,276 bytes, owner's download 2026-09-09). **The
file name is wrong: the PDF is Chan & Reineccius 1994, the RSC book chapter (pp. 131-139), not anything
from 2005.** The name is left as it is; every reference to it in the repository should be to
"Chan 1994b" / this dossier. The PDF is an Elsevier-produced scan (metadata: "Acrobat 4.0 Import
Plug-in", subject "Maillard Reaction in Chemistry, Food, and Health (xx) 131-139") with an OCR text
layer; read from `scratchpad/articles/chan2005.txt`, with pages 136-137 (Figures 1-3) rendered at 110
dpi to check the axes. OCR quality: digits in the two tables and the loadings are clean and the mol-g
pairs check against the molar masses (section 2); Greek letters, the degree sign and some ligatures are
mangled ("11YC" for 115 C, "kcallmole"). Tables I and II re-typed in full below. Figures 1-3 are bar/line
charts on a "Conc. in reaction mixture (mole)" ordinate and are FIGURE-ONLY. There is no supplementary
material for a 1994 book chapter. The companion chapter on the sulfur volatiles from the same run is on
disk as `chan1994_extraction.md` (ACS Symp. Ser. 564, DOI 10.1021/bk-1994-0564.ch010); that dossier's
audit of the reactor, buffer and quantification applies here verbatim, because the data come from the same
five-temperature, three-pH run.

## 0. Identity

| field | value |
|---|---|
| Title | "The Reaction Kinetics for the Formation of Isovaleraldehyde, 2-Acetyl-1-pyrroline, di(H)di(OH)-6-Methylpyranone, Phenylacetaldehyde, 5-Methyl-2-phenyl-2-hexenal, and 2-Acetylfuran in Model Systems" |
| Authors | F. Chan and G. A. Reineccius, Department of Food Science and Nutrition, University of Minnesota, 1334 Eckles Avenue, St. Paul, MN 55108 |
| Book | *Maillard Reactions in Chemistry, Food and Health*, T. P. Labuza, G. A. Reineccius, V. M. Monnier, J. O'Brien, J. W. Baynes (eds.), The Royal Society of Chemistry, Cambridge, 1994 (series number not printed in the PDF) |
| Pages | 131-139 (p. 131 title and summary; Table I p. 134; Table II p. 135; Figures 1-3 pp. 136-137; conclusions and CAS numbers p. 137; references pp. 137-139) |
| DOI | **none printed anywhere in the PDF.** Do not attach one from memory; the repository citation is book + pages |
| PDF file name | `data/articles/chan2005.pdf` — **misnamed** (see above); the pdfinfo Title field carries the correct chapter title |
| On disk vs SI | the whole chapter is on disk; no supplementary material exists |
| Naming | isovaleraldehyde = 3-methylbutanal (registry `3_methylbutanal`, alias present); phenylacetaldehyde = "2-phenylethanal" in the figure legends; di(H)di(OH)-6-methylpyranone = DDMP = 2,3-dihydro-3,5-dihydroxy-6-methyl-4H-pyran-4-one; "Chan and Reineccius 1994" in this chapter's own reference list is the sulfur chapter (`chan1994_extraction.md`), then "in press" |
| CAS printed by the authors (p. 137) | isovaleraldehyde 590-86-3; 2-acetyl-1-pyrroline 85213-22-5; phenylacetaldehyde 122-78-1; 4-methylthiazole 693-95-8; 5-methyl-2-phenyl-2-hexenal 21834-92-4; 2-acetylfuran 1192-62-7 (all correct; DDMP has none) |
| Who quotes it | Cremer & Eichner 2000 ("80.4 kJ/mol" for 3-MB in water; `cremer2000_extraction.md`), Balagiannis 2009 (same 80.4; `balagiannis2009_extraction.md`), Huang 2017, Parker 2013 (ref 16; `parker2013_extraction.md` Flag 2 warned the chapter was not on disk), Balagiannis 2015 ("60 to 129 kJ/mol" for both 1994 chapters together) |

## 1. Why it matters

The amino-acid-identity wave planned in `results/validation/kinetic_core_b19_prereg_draft.md` needs, per
amino acid, a Strecker-aldehyde rate at two or more temperatures in water; its row table names this
chapter as "the OTHER 1994 chapter: 3-methylbutanal and phenylacetaldehyde pseudo-zero-order rates and
barriers at pH 6 to 8, 75 to 115 C". Half of that is true. The chapter prints **barriers** for the leucine
and phenylalanine aldehydes (and for proline's 2-acetyl-1-pyrroline, the ring the draft lists under "to
find"), all from one pot that contains the four amino acids at once, in the temperature window the trunk
runs in (75-115 C against the trunk's 80-120 C). It prints **no rate constants** — not in a table, not
in a figure, not in the text — so nothing here fixes the magnitude of the Strecker step `k_strecker` (rule
R07) for any amino acid; the pyrazine step's two constants (`parameters_pyrazine.py`, `FROZEN_B18`:
log10 k = -6.54 and -7.53 L/(mmol min) at 100 C, Ea 103.1 and 114.9 kJ/mol, second order in dicarbonyl
and glycine) remain the only measured Strecker magnitudes on the trunk. What the chapter does give the
wave is (a) a within-study barrier ordering across three amino acids in the same pot — proline's product
60 kJ/mol < leucine's 80 kJ/mol < phenylalanine's 90 kJ/mol — which is a test of the draft's assumption
that one pH term and one barrier can be shared across amino acids, and (b) the observation that
phenylacetaldehyde and isovaleraldehyde are consumed by an aldol condensation at pH 7-8 fast enough to
be visible at room temperature, which means any per-amino-acid Strecker rate fitted on a net aldehyde
level in a multi-amino-acid pot is a net-of-sink rate.

## 2. Methods as they matter to a model

- **Pot.** "Glucose (0.5 mole = 90 g), methionine (0.075 mole = 11.19 g), phenylalanine (0.075 mole =
  12.39 g), proline (0.075 mole = 8.64 g) and leucine (0.075 mole = 9.84 g) were dissolved in 400 mL of
  distilled water and pH adjusted to 6, 7 or 8 with a 0.1M phosphate buffer and the appropriate amount of
  NaOH" (the OCR reads "6.7 or 8"; the results discuss pH 6, 7 and 8 throughout, and the summary's "pH
  6.7 & 8" is the same three values). All five mol-to-gram pairs check against the molar masses (glucose
  180.16, Met 149.21, Phe 165.19, Pro 115.13, Leu 131.17 g/mol), so the loadings are absolute moles.
  **Molarities (mine): glucose 0.5 mol / 0.400 L = 1250 mmol/L; each amino acid 187.5 mmol/L; total amino
  acid 750 mmol/L; glucose : total amine 1.67 : 1; glucose : leucine 6.67 : 1.** Buffer: "0.1M
  phosphate" is the concentration of the stock used to adjust pH, not a stated final molarity; whatever
  was added, 0.1 mol/L phosphate against 750 mmol/L of titratable amine and 1.25 mol/L of sugar does not
  hold pH (the sulfur dossier's reading: **pH 6 / 7 / 8 are initial values of a drifting system**). The
  pH of every sample was measured after cooling and none is printed.
- **Reactor.** 600 mL Parr 4563 pressure reactor with a Parr 4842 controller, filled with the entire
  400 mL charge, sealed, stirred; "initial temperature was noted and start time determined when the
  solutions reached reaction temperatures" (heat-up time not printed). 50 mL samples withdrawn through
  the sampling port on the schedule below; about 50 mL of ultra-high-purity N2 added after each sampling
  to restore pressure. Over a run 250 mL of the 400 mL charge is withdrawn, so the headspace grows from
  about 200 to about 450 mL of N2 along each time series (the sulfur dossier's confounder; for the
  aldehydes here it changes the liquid-headspace partition of isovaleraldehyde, b.p. 92 C, along the
  fitted axis).
- **Schedule (printed, p. 132):** 75 C: 1.5, 3.0, 4.5, 6.0, 7.5 h; 85 C: 0.5, 1.0, 1.5, 2.0, 2.5 h; 95 C:
  20, 40, 60, 100, 120 min; 105 C: 10, 20, 30, 40, 50 min; 115 C: 5, 10, 15, 20, 25 min. Five points per
  isotherm, five temperatures, three pH, in duplicate ("Duplicates at each temperature and pH series
  was done"). A zero-time extract exists only for the pH 7 series ("to serve as an anchor point").
- **Work-up.** Sample cooled to room temperature in an ice bath, pH measured, extracted 3 x 5 mL
  dichloromethane containing **500 ppm 4-methylthiazole as internal standard**, dried over MgSO4,
  filtered, concentrated to 0.5-1.0 mL under N2.
- **Identification.** HP 5890 / 5970 MSD, 70 eV, DB-5 30 m x 0.32 mm x 1 µm, 40 C (3 min) to 250 C at
  5 C/min, 20 : 1 split, 1 µL; NBS library + published retention indices + co-chromatography with
  authentic compounds.
- **Quantification ("to collect kinetic data").** HP 5890 Series II with an HP 5921A atomic emission
  detector, carbon 193 nm, nitrogen 174 nm, sulfur 181 nm channels, same column, 45 : 1 split, 1 µL.
  **No calibration curve, no response factor, no recovery, no LOD is printed.** The kinetic quantity is
  therefore an internal-standard-normalised AED channel response; the sulfur chapter labels its figures
  "Amount/ppm" on an undefined basis, and this chapter's only printed level is "151 ppm" (below). Which
  channel served the non-sulfur, non-nitrogen compounds (isovaleraldehyde, phenylacetaldehyde,
  2-acetylfuran, DDMP, the hexenal: carbon only) against a sulfur-and-nitrogen internal standard is not
  stated. **Absolute levels are not transportable; ratios within a compound across temperature are.**
- **Kinetic analysis.** Labuza's "Water Analyzer Series - Reaction Kinetics Program Version 2.09";
  order chosen by the r2 of quantity-versus-time regressions ("the r2 for all the compounds except
  5-methyl-2-phenyl-2-hexenal at zero order were consistently better than any other order"); k = slope;
  Arrhenius on the five temperatures; **Table I is "average" Ea** — the average over pH (6, 7, 8) and
  over the two replicates, six numbers from what were 36 ladders. No per-pH, per-replicate Ea, no r2, no
  confidence interval is printed here (contrast the sulfur chapter, which prints replicate-by-replicate Ea
  and r2). Analysis of covariance (MacAnova 3.1): pH significant for the rate of all compounds except
  2-acetylfuran, phenylacetaldehyde and the hexenal; temperature x pH not significant; replicate effect
  significant for 2-acetylfuran and the hexenal.
- **Second experiment (Figures 1-3).** "Reaction of isovaleraldehyde and 2-phenylethanal at 95 C" at pH
  6, 7, 8, sampled at 20, 40, 60 min (pH 7 also at 0): the two aldehydes heated together, presumably in
  the same buffer, and the aldol product 5-methyl-2-phenyl-2-hexenal followed. **The set-up of this
  experiment (amounts charged, volume, whether glucose or amino acids were present) is not described
  anywhere in the chapter.** The ordinate is "Conc. in reaction mixture (mole)", of order 1e-3.
- **Unit conversions used below.** 1 kcal/mol = 4.184 kJ/mol.

## 3. Tables re-typed

### Table I. "Average activation energies for the Strecker aldehydes and Maillard reaction products" (p. 134)

| Compound | Activation Energy (kcal/mole), as printed | kJ/mol (mine, x 4.184) |
|---|---:|---:|
| 2-Acetyl-1-pyrroline | 14.4 | 60.2 |
| di(H)di(OH)-6-Methylpyranone | 16.1 | 67.4 |
| 2-Acetylfuran | 17.7 | 74.1 |
| Isovaleraldehyde | 19.2 | 80.3 |
| Phenylacetaldehyde | 21.5 | 90.0 |
| 5-Methyl-2-phenyl-2-hexenal | 23.1 | 96.7 |

No footnote, no SD, no n, no r2, no unit for k, no per-pH breakdown. The hexenal's Ea is on a
**first-order** fit (Conclusions); the other five on zero order.

### Table II. "Comparison of our activation energy data with Schirle-Keller and Reineccius (1990)" (p. 135)

| | Published (1) | Current study |
|---|---:|---:|
| 2-Acetylfuran | 36.2 kcal/mole | 17.7 kcal/mole |
| di(H)di(OH)-6-Methylpyranone | 30.7 kcal/mole | 16.0 kcal/mole |

(1) data of Schirle-Keller and Reineccius (1990) [ACS Symp. Ser. 490, glucose + cysteine, unbuffered].
Note DDMP is 16.0 here and 16.1 in Table I (Flag 3). In kJ/mol: 151.5 / 74.1 and 128.4 / 66.9.

### Numbers in the running text (everything else is figure-only)

- Summary: "activation energies ranged from 14.3 to 23.1 kcal/mole" (59.8-96.7 kJ/mol); Table I's
  minimum is 14.4; the text (p. 134) and Conclusions say "14 to 23".
- Comparators quoted by the authors: Schirle-Keller & Reineccius 1990, oxygen heterocycles "28 to 33
  kcal/mole" (117-138 kJ/mol); Leahy & Reineccius 1989a,b, pyrazines "27 to 45 kcal/mole" (113-188
  kJ/mol); Leahy 1985 thesis: pyrazine Ea highest at pH 5, lowest at pH 7 (three temperatures, three
  points each, "inadequate for statistical calculations").
- **"a zero time extract had been taken for the pH 7 sample ... at room temperature and zero time,
  already 151 ppm of 5-methyl-2-phenyl-2-hexenal was found"** — the one printed level in the chapter;
  ppm basis undefined (Flag 6). It refers to the second experiment (Figure 2, the pH 7 panel, t = 0).
- "the concentration of 5-methyl-2-phenyl-2-hexenal reached a maximum at 20 min heating and decreased
  with continued heating" (second experiment, 95 C).
- Qualitative pH statements from the main run: isovaleraldehyde higher at pH 7 than pH 8;
  phenylacetaldehyde highest at pH 6; the hexenal higher at pH 7 and 8 than at 6; 2-acetyl-1-pyrroline
  rising with time and pH with "little difference between pHs 6 and 7"; 2-acetylfuran rising with
  temperature, pH and time; DDMP "only began to occur at temperatures of 95 C or higher".
- Second experiment: "the reaction/degradation of phenylacetaldehyde was lowest at pH 6 and increased
  with increasing pH"; "the concentration of isovaleraldehyde in the reaction mixture did not vary with
  pH's" (Figures 1-3).

### Figures 1-3 (pp. 136-137) — FIGURE-ONLY

Three panels (pH 6, 7, 8) at 95 C: bars for 5-methyl-2-phenyl-2-hexenal, lines for isovaleraldehyde and
2-phenylethanal, at 20 / 40 / 60 min (pH 7 also 0 min); ordinate "Conc. in reaction mixture (mole)",
scale 0 to 3e-3 or 4e-3. No number is typed from them. Shape: at pH 7 and 8 phenylacetaldehyde is near
its floor by 20 min while isovaleraldehyde stays high; at pH 6 both stay high and the hexenal is low.

## 4. Kinetic numbers the repository can use

Registry mapping: isovaleraldehyde -> `3_methylbutanal` (alias "isovaleraldehyde" present);
phenylacetaldehyde -> `phenylacetaldehyde`; 2-acetylfuran -> `2_acetylfuran`; **2-acetyl-1-pyrroline,
DDMP, 5-methyl-2-phenyl-2-hexenal, 4-methylthiazole (the internal standard; only its 2-alkyl and
dihydro relatives are keyed), glucose, leucine, phenylalanine, proline, methionine -> not in
`data/keys/compounds.yml`.**

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| Strecker of leucine (R07 + supply), net of sinks | apparent Ea, isovaleraldehyde formation | 19.2 (= 80.3 kJ/mol) | kcal/mol | glucose 1250 + Met, Phe, Pro, Leu 187.5 each mmol/L, phosphate (nominal 0.1 mol/L) pH 6-8 initial, 75-115 C, N2 headspace, stirred Parr reactor; average over pH 6/7/8 and two replicates | pseudo-zero order in time (k = slope of AED response vs t) | Table I, p. 134 | measured_barrier (average; no SD; response-normalised, so unit-free) |
| Strecker of phenylalanine, net of the aldol sink | apparent Ea, phenylacetaldehyde formation | 21.5 (= 90.0 kJ/mol) | kcal/mol | same | pseudo-zero order | Table I | measured_barrier (average; no SD) |
| proline route to 2-acetyl-1-pyrroline | apparent Ea | 14.4 (= 60.2 kJ/mol) | kcal/mol | same (multi-amino-acid pot; the authors say secondary products are "very system dependent") | pseudo-zero order | Table I | measured_barrier (average; no SD) |
| sugar path, 2-acetylfuran | apparent Ea | 17.7 (= 74.1 kJ/mol) | kcal/mol | same | pseudo-zero order | Tables I, II | measured_barrier |
| sugar path, DDMP | apparent Ea | 16.1 (Table I) / 16.0 (Table II) (= 67.4 / 66.9 kJ/mol) | kcal/mol | same; DDMP appears only at >= 95 C, so the ladder is at most three temperatures | pseudo-zero order | Tables I, II | measured_barrier (weak: see conditions) |
| aldol, isovaleraldehyde + phenylacetaldehyde -> hexenal | apparent Ea | 23.1 (= 96.7 kJ/mol) | kcal/mol | same | first order | Table I | measured_barrier (on a first-order fit of a product that peaks at 20 min at 95 C; meaning unclear) |
| within-study barrier ordering | Ea(2-AP) : Ea(3-MB) : Ea(PAA) | 60.2 : 80.3 : 90.0 | kJ/mol | one pot, one detector, one analysis | — | Table I | within_study_ratio (of barriers, not rates) |
| any rate constant, any compound, any T or pH | — | **not printed** | — | — | — | whole chapter | — |
| hexenal at t = 0, pH 7 second experiment | 151 | ppm (basis undefined) | room temperature, before heating | — | text p. 135 | level_only (peak_area_only in substance: AED response vs a 500 ppm IS, no response factor) |
| aldehyde and hexenal levels vs time at 95 C, pH 6/7/8 | — | "mole" | second experiment, set-up undescribed | — | Figures 1-3 | figure_only |
| pH dependence of the six rates | significant for all but 2-acetylfuran, PAA, hexenal (ANCOVA) | — | pH 6-8 initial | — | text p. 134 | level_only (a significance statement, no effect size) |

**Can a second-order constant in water be derived?** No. Nothing to convert: no k is printed. Had the
zero-order k been printed (as an AED response per minute) they would still hide glucose 1250 mmol/L,
leucine 187.5 mmol/L and 750 mmol/L of competing amine, and lack a response factor, so even then only
the ratio form k(T1)/k(T2) — the barrier — would transport. **What transports: the barriers, as
whole-cascade apparent barriers from a sugar + amino-acid pot in water, initial pH 6-8, 75-115 C.**
They are not the barrier of the dicarbonyl + amino-acid step alone; the authors say so ("the kinetics
observed are the dynamic result of formation vs consumption").

**Comparison with the trunk's constants.** The pyrazine step's measured Strecker barriers on fed
dicarbonyls are 103.1 kJ/mol (glyoxal + alanine) and 114.9 kJ/mol (methylglyoxal + alanine)
(`FROZEN_B18`); Martins' supply steps are 125.0 +/- 4.7 (Amadori -> methylglyoxal) and 107 +/- 7.3
(Amadori -> 1-deoxyglucosone) kJ/mol (`MARTINS_M4`); Cremer & Eichner's low-moisture barriers are
115-124 kJ/mol; Balagiannis 2009's glucose-to-intermediate barrier in liver extract is 137 +/- 15
kJ/mol. Chan's 80 (leucine) and 90 (phenylalanine) kJ/mol sit 20-40 kJ/mol below all of these. Three
readings are possible and the chapter cannot separate them: (i) the response-normalised ppm scale is
not linear in concentration for the carbon channel, flattening the ladder; (ii) the aldehyde sinks
(aldol, further Maillard reactions) grow faster with temperature than the source, lowering the net
apparent barrier; (iii) at 1.25 mol/L glucose and 0.75 mol/L amine the supply of dicarbonyl is not the
limiting step it is in a 200 + 200 mmol/L pot. Reading (ii) is supported inside the paper by the
hexenal experiment. Record 80.3 kJ/mol as a lower bound for the leucine aldehyde's apparent barrier
in water, not as the step's barrier.

## 5. Flags

1. **File name.** `chan2005.pdf` is Chan & Reineccius 1994 (RSC, pp. 131-139). Not renamed; cite as
   Chan 1994b. The on-disk `chan1994_extraction.md` is the ACS sulfur chapter and must not be cited for
   isovaleraldehyde or phenylacetaldehyde.
2. **No rate constants anywhere.** Table I has Ea only; there is no k table, no k figure, no
   Arrhenius plot. The 36 ladders (6 compounds x 3 pH x 2 replicates) are reduced to six averaged
   numbers without dispersion. Request from the authors (if the 1994 data survive): the k tables per
   compound / pH / replicate in the AED response unit, and the per-pH Ea with 95 % limits that the text
   says were computed ("when the 95 % confidence limits are taken into account").
3. **Internal inconsistencies of the printed numbers.** Summary "14.3 to 23.1" vs Table I minimum 14.4;
   DDMP 16.1 (Table I) vs 16.0 (Table II); Summary says "pH 6.7 & 8" and "75 to 11YC" (OCR of 115 C).
   The 80.4 kJ/mol quoted by Cremer 2000, Balagiannis 2009 and Huang 2017 is 19.2 x 4.184 = 80.33
   rounded up; it is the same number.
4. **The pot is a four-amino-acid competition.** Isovaleraldehyde and phenylacetaldehyde form beside
   methionine and proline products at 187.5 mmol/L each; the authors themselves warn that secondary
   products (2-acetyl-1-pyrroline) and sugar fragments (2-acetylfuran, DDMP) are "very system
   dependent" in such a pot. The Strecker aldehyde barriers are "likely less" affected, by the authors'
   judgement, not by a control.
5. **pH is an initial value.** "0.1M phosphate" against 750 mmol/L amine; the measured sample pHs are
   not printed. The "pH 6 / 7 / 8" columns of the ANCOVA are labels for drifting series. Any pH term
   read from this chapter (the ordering PAA: pH 6 highest; 3-MB: pH 7 > 8; 2-AP: 8 > 7 ~ 6) is a
   statement about initial pH.
6. **Quantification is AED response against a 500 ppm sulfur-nitrogen internal standard with no
   response factors.** The "151 ppm" of hexenal at t = 0 and every figure ordinate are on an undefined
   ppm / "mole" basis. Nothing here is a concentration. Class everything level-like as
   peak_area_only in substance.
7. **Figures 1-3 belong to an undescribed experiment.** The methods describe only the glucose +
   four-amino-acid run; the figures follow isovaleraldehyde + phenylacetaldehyde heated together at 95
   C. Charged amounts, volume, buffer and the presence or absence of glucose are not stated. The
   ordinate "Conc. in reaction mixture (mole)" reads as an absolute amount in the reactor, not a
   molarity.
8. **The zero-order claim is a five-point r2 comparison** on a quantity that is itself net of sinks;
   the hexenal peaks at 20 min at 95 C, so at least one "zero-order" product is in its consumption
   phase inside the fitted window. Five points, five temperatures, two replicates; the sulfur chapter's
   audit of the same run found leg-to-leg Arrhenius slope breaks >= 2x in every ladder.
9. **Headspace grows along each series** (250 of 400 mL withdrawn; N2 top-ups of ~50 mL). For
   isovaleraldehyde (b.p. 92 C, Henry constant unfavourable to the liquid at 95-115 C) the liquid
   concentration at late points is biased low; a zero-order slope is biased low and the barrier with it.
10. **Temperature window.** The 75 C series ran to 7.5 h and the 115 C series to 25 min; heat-up time of
    400 mL in a Parr reactor to 115 C is not printed and "start time determined when the solutions
    reached reaction temperatures" means the early points at 105-115 C carry a warm-up contribution.
11. **Registry gaps** (`data/keys/compounds.yml`): 2-acetyl-1-pyrroline, DDMP,
    5-methyl-2-phenyl-2-hexenal, 4-methylthiazole and all five reactants are unkeyed. The planned
    proline route needs a `2_acetyl_1_pyrroline` key before any row can be written.
12. **What to fetch next for the same question.** Schirle-Keller & Reineccius 1990/1992 (ACS 490) for
    the 2-acetylfuran / DDMP comparison; Stahl & Parliment 1994 (ACS 543, pp. 251-262) for
    proline-glucose high-temperature short-time kinetics (2-acetyl-1-pyrroline); Leahy 1985 thesis
    for the per-pH pyrazine ladders behind `PYRAZINE_PH_SLOPES`.
