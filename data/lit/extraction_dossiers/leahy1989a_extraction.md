# Leahy & Reineccius 1989 (ch. 18) — EXTRACTION (0.1 M lysine + 0.1 M glucose at pH 5.0 / 7.0 / 9.0, 75/85/95 C, up to 24 h; and nonfat dry milk at aw 0.32-0.84, 95 C, 0.5-3 h)
### Part Two of the pair: the pH ladder for pyrazine, methylpyrazine and dimethylpyrazine rates, and a water-activity series in a milk powder.

**Source on disk:** `data/articles/leahy1989a.pdf` (13 pp., owner's download, 2026-09-08). ACS
Symposium Series 409 (*Thermal Generation of Aromas*, Parliment, McGorrin, Ho eds.), chapter 18,
pp. 196-208. The text layer is upper-case OCR from p. 196 onward but numerically clean; **Table I
was verified cell by cell against a 130-dpi raster of printed p. 200** (`scratchpad/img/l89a-05.png`);
Tables II-V were taken from the text layer (unambiguous). Figures 1-4 are FIGURE-ONLY and were not
read.

## 0. Identity

| field | value |
|---|---|
| Title | "Kinetics of the Formation of Alkylpyrazines. Effect of pH and Water Activity" |
| Authors | M. M. Leahy (Ocean Spray Cranberries, Lakeville-Middleboro MA) and Gary A. Reineccius (Univ. of Minnesota, St. Paul) |
| Venue | ACS Symposium Series 409, *Thermal Generation of Aromas*, 1989, ch. 18, pp. 196-208. Received May 31 1989; published October 3 1989 |
| DOI | 10.1021/bk-1989-0409.ch018 |
| Naming | "NFDM" = nonfat dry milk (OCR renders it "NFEM"); "2-methylpyrazine" = methylpyrazine; "Part One" = ch. 7 of ACS Symp. Ser. 388 (`leahy1989_extraction.md`; the chapter's ref 18 cites it with the wrong page 76 as "p. 76") |
| Companion | The pH 9.0 rows of Table I and the pH 9.0 column of Table III are the lysine-glucose data of ch. 7, re-printed (with small discrepancies, §5). What this chapter ADDS: pH 7.0 and pH 5.0 arms at the same three temperatures; Ea per pH; the 2 h / 95 C yield ladder 0.0239 / 6.19 / 13.1 ppm; the NFDM aw series |
| Cited in repo | `apriyantono1993_extraction.md` uses Table III's pH 5.0 column (pyrazine 63.2 %, methylpyrazine 36.8 %, quoted there as 632 / 368 g/kg) as a hexose comparator; `chan1994_extraction.md` items 74-75 |

## 1. Why it matters

The pyrazine lane needs a pH law. This is the only source in the corpus with **the same
sugar-amino acid system, method and temperatures at three pH values two units apart, with printed
rates and an Ea at each pH**. It gives the within-study ratios k(pH 9)/k(pH 7) ~ 2-3 and
k(pH 9)/k(pH 5) ~ 20-60 for pyrazine, and the observation that the dimethylpyrazines fall below
quantification at pH 5. The aw part gives the shape of a rate-vs-aw curve (maximum near aw 0.75)
in a real matrix, but in units that are not concentrations. Limits as for ch. 7: 75-95 C, lysine
only, borate/citrate-phosphate buffers, conversions of order 0.1 %, 3-point Arrhenius.

## 2. Methods as they matter to a model

- **pH study reactants.** "ten ml of buffered solutions at pH 5.0, 7.0 and 9.0 of 0.1M L-lysine
  monohydrochloride and 0.1M D-glucose ... were heated in Teflon-capped 25 mm (o.d.) x 150 mm Pyrex
  test tubes in a water bath at 75, 85, and 95 C for up to 24 hr." So **100 mmol/L lysine·HCl +
  100 mmol/L glucose**, 100 mmol/L chloride present.
- **Buffers.** "Citrate-phosphate buffers (0.1M) were used to achieve a pH of 5.0 and 7.0 and borate
  buffer (0.1M) was used for a pH of 9.0." The buffer changes between pH 7 and 9; the authors
  assert "a change in buffer type did not affect pyrazine formation" without printing a test.
- **Sampling.** "Samples were taken at 7 to 8 time intervals. Eighteen to 22 total samples per
  temperature and pH were analyzed." Two zero points per regression; duplicates early, triplicates
  late. Sampling times not printed.
- **Work-up and analysis.** Identical to ch. 7: after heating, adjust to pH 9.0 with 0.1 N NaOH; add
  1 mL of 2 ppm 2-methoxypyrazine (IS); final volume 15 mL; HP 7675A purge-and-trap onto Tenax;
  HP 5880A GC with nitrogen-phosphorus detector; "Specifics of sample preparation and chromatographic
  analysis have been described previously (18)". Quantification by the IS method with empirical
  response factors (ch. 7 §2); unit ppm = µg/mL; 10 -> 15 mL dilution basis unstated (ch. 7 Flag 3).
  Identification by co-chromatography with standards and GC-MS.
- **aw study.** Fresh NFDM (Maple Island), five 30 g samples equilibrated 2 weeks at room temperature
  in evacuated desiccators over saturated MgCl2, NaBr, NaCl, KCl (nominal aw 0.32, 0.58, 0.75, 0.85);
  measured initial aw "0.319, 0.583, 0.747 and 0.841" (Kaymont-Rotronics hygrometer). Triplicate 5 g
  samples in 75 mm o.d. x 150 mm Pyrex culture tubes with Teflon-sealed caps, 95 C water bath,
  "five sampling times, ranging from 0.5 to 3 hr". Powder broken up and **purged directly as a dry
  powder**; no internal standard possible; quantification "relative to response factors determined
  for external standard solutions"; then "standardized concentration units" = relative
  concentration x 100 / initial dry weight of sample (moisture by the GC method of Reineccius &
  Addis 1973). **The Table IV "ppm/hr" are therefore not solution concentrations** (Flag 2).
- **Kinetic treatment (verbatim).** "dA/dt = k A^n ... A = A0 + kt ... The formation of pyrazines fit
  a zero order reaction. Plotting concentrations of pyrazines formed versus time of reaction gave the
  better fit of the line, usually with a coefficient of determination (r2) of greater than 0.95."
  A in ppm, t in h, k in ppm/h; Arrhenius with R = 1.986 cal/mol/K, Ea in kcal/mol, three
  temperatures. "Since rate constants were determined at only 3 temperatures, only 3 data points were
  used to determine activation energies. Data at other temperatures is necessary to make any further
  comparisons among activation energies."

## 3. Tables re-typed

### Table I. "Regressions for the formation of pyrazines (0.1M lysine-glucose systems)"

Columns: compound / pH, temperature, k (ppm/hr), intercept, n (number of data points), r2. Added:
k in µmol L-1 min-1 (= ppm/h / MW x 1000 / 60; MW pyrazine 80.088, methylpyrazine 94.115,
dimethylpyrazines 108.141), on the assumption that ppm is the reaction-mixture concentration.

| compound | pH | T (C) | k (ppm/h) | intercept (ppm) | n | r2 | k (µmol L-1 min-1) |
|---|---:|---:|---:|---:|---:|---:|---:|
| pyrazine | 9.0 | 95 | 3.596 | 0.0596 | 22 | 0.994 | 0.748 |
| | 9.0 | 85 | 0.490 | 0.458 | 22 | 0.960 | 0.102 |
| | 9.0 | 75 | 0.214 | 0.279 | 20 | 0.965 | 0.0445 |
| | 7.0 | 95 | 1.346 | 0.762 | 19 | 0.962 | 0.280 |
| | 7.0 | 85 | 0.159 | 0.609 | 22 | 0.899 | 0.0331 |
| | 7.0 | 75 | 0.0957 | -0.0014 | 22 | 0.989 | 0.0199 |
| | 5.0 | 95 | 0.0938 | -0.0647 | 22 | 0.984 | 0.0195 |
| | 5.0 | 85 | 0.0232 | -0.0586 | 22 | 0.966 | 0.00483 |
| | 5.0 | 75 | 0.00356 | -0.0185 | 17 | 0.928 | 0.000741 |
| 2-methylpyrazine | 9.0 | 95 | 2.837 | -0.104 | 22 | 0.995 | 0.502 |
| | 9.0 | 85 | 0.422 | 0.142 | 22 | 0.967 | 0.0747 |
| | 9.0 | 75 | 0.159 | 0.091 | 20 | 0.941 | 0.0282 |
| | 7.0 | 95 | 1.367 | -0.070 | 19 | 0.981 | 0.242 |
| | 7.0 | 85 | 0.276 | 0.204 | 22 | 0.967 | 0.0489 |
| | 7.0 | 75 | 0.0945 | 0.028 | 22 | 0.981 | 0.0167 |
| | 5.0 | 95 | 0.00636 | -0.00212 | 22 | 0.890 | 0.00113 |
| | 5.0 | 85 | 0.00279 | -0.00430 | 22 | 0.912 | 0.000494 |
| 2,5-dimethylpyrazine | 9.0 | 95 | 0.186 | -0.0457 | 20 | 0.976 | 0.0287 |
| | 9.0 | 85 | 0.0247 | -0.00604 | 22 | 0.985 | 0.00381 |
| | 9.0 | 75 | 0.00668 | -0.00536 | 16 | 0.942 | 0.00103 |
| | 7.0 | 95 | 0.0630 | -0.0051 | 19 | 0.965 | 0.00971 |
| | 7.0 | 85 | 0.00949 | -0.00229 | 22 | 0.974 | 0.00146 |
| | 7.0 | 75 | 0.00209 | -0.00574 | 14 | 0.857 | 0.000322 |
| 2,3-dimethylpyrazine | 9.0 | 95 | 0.0229 | -0.0057 | 16 | 0.948 | 0.00353 |
| | 9.0 | 85 | 0.00309 | -0.00017 | 18 | 0.978 | 0.000476 |
| | 9.0 | 75 | 0.000677 | -0.00086 | 12 | 0.958 | 0.000104 |

No rows exist for 2-methylpyrazine at pH 5.0 / 75 C, for either dimethylpyrazine at pH 5.0, or for
2,3-dimethylpyrazine at pH 7.0: "quantification of the dimethylpyrazines could only be accomplished at
pH's greater than 5.0". 2,6-dimethylpyrazine does not appear (not detected in lysine-glucose, per ch. 7).

**Within-study pH ratios (mine, from Table I):**

| compound | ratio | 95 C | 85 C | 75 C |
|---|---|---:|---:|---:|
| pyrazine | k(pH 9)/k(pH 7) | 2.67 | 3.08 | 2.24 |
| pyrazine | k(pH 9)/k(pH 5) | 38.3 | 21.1 | 60.1 |
| pyrazine | k(pH 7)/k(pH 5) | 14.4 | 6.9 | 26.9 |
| 2-methylpyrazine | k(pH 9)/k(pH 7) | 2.08 | 1.53 | 1.68 |
| 2-methylpyrazine | k(pH 9)/k(pH 5) | 446 | 151 | — |
| 2-methylpyrazine | k(pH 7)/k(pH 5) | 215 | 99 | — |
| 2,5-dimethylpyrazine | k(pH 9)/k(pH 7) | 2.95 | 2.60 | 3.20 |

Text: regressing k on pH "gave a good fit of the line, with R2 values of 0.974 and 0.999" for
pyrazine and 2-methylpyrazine (Figure 1; which temperature's k, and the slope, are not printed —
FIGURE-ONLY).

### Table II. "Activation energies for formation of pyrazines (0.1M lysine-glucose systems)" (Ea in kcal/mol)

kJ/mol (x 4.184) and my 3-point refit of Table I added.

| compound | pH | Ea printed (kcal/mol) | Ea (kJ/mol) | refit from Table I (kJ/mol) |
|---|---:|---:|---:|---:|
| pyrazine | 9.0 | 35.8 | 149.8 | 149.7 |
| | 7.0 | 33.4 | 139.7 | 140.0 |
| | 5.0 | 41.8 | 174.9 | 174.5 |
| 2-methylpyrazine | 9.0 | 36.7 | 153.6 | 153.0 (ch. 7 prints 36.6) |
| | 7.0 | 34.0 | 142.3 | 142.1 |
| 2,5-dimethylpyrazine | 9.0 | 41.9 | 175.3 | 176.8 (ch. 7 prints 42.3) |
| | 7.0 | 43.7 | 182.8 | 181.2 |
| 2,3-dimethylpyrazine | 9.0 | 44.8 | 187.4 | 187.3 |

Text: "Activation energies for alkylpyrazine formation ... ranging from 33 to 45 kcal/mole ...
Activation energies for pyrazine and 2-methylpyrazine formation were approximately 35 kcal/mole,
while those for the dimethylpyrazines were slightly higher." "For pyrazine, the activation energy is
lowest at pH 7.0 at 33 kcal/mole." (The pH 5.0 2-methylpyrazine pair, 2 points, would give 90
kJ/mol — not printed by the authors and not to be used.)

### Table III. "Effect of pH on distribution pattern of pyrazines — 0.1M lysine-glucose, 2 hr, 95 C"

| compound | pH 5.0 | pH 7.0 | pH 9.0 |
|---|---:|---:|---:|
| pyrazine (%) | 63.2 | 58.7 | 55.8 |
| 2-methylpyrazine (%) | 36.8 | 39.6 | 41.5 |
| 2,5-dimethylpyrazine (%) | — | 1.7 | 2.4 |
| 2,3-dimethylpyrazine (%) | — | — | 0.3 |
| **TOTAL (ppm)** | **0.0239** | **6.19** | **13.1** |

Column sums 100.0 / 100.0 / 100.0. Text: "a total of 13,000 ppb pyrazines produced at pH 9.0 and
only 24 ppb at pH 5.0". Yield ratios: pH 9 / pH 7 = 2.12; pH 9 / pH 5 = 548; pH 7 / pH 5 = 259.
Molar total at pH 9 = 0.152 mmol/L (ch. 7 §3 arithmetic), i.e. 0.3 % of the lysine.

### Table IV. "Regressions for the formation of pyrazines in NFDM as a function of time (ppm/hr) at 95 C at aws 0.32, 0.58, 0.75 and 0.84"

y = "standardized concentration" (relative concentration x 100 / initial dry weight; see §2), x =
time (h). n = number of data points.

| compound | aw (measured) | regression | n | R2 |
|---|---:|---|---:|---:|
| pyrazine | 0.32 (0.319) | y = 0.324x + 0.00856 | 17 | 0.895 |
| | 0.58 (0.583) | y = 0.681x - 0.110 | 17 | 0.930 |
| | 0.75 (0.747) | y = 1.371x - 0.102 | 17 | 0.971 |
| | 0.84 (0.841) | y = 1.318x + 0.0307 | 18 | 0.940 |
| 2-methylpyrazine | 0.32 | y = 0.271x - 0.0199 | 17 | 0.922 |
| | 0.58 | y = 0.626x - 0.135 | 17 | 0.899 |
| | 0.75 | y = 0.886x - 0.133 | 17 | 0.975 |
| | 0.84 | y = 0.718x - 0.0838 | 18 | 0.955 |

Only pyrazine and 2-methylpyrazine were identified in NFDM. "Reaction rates increased with aw over
the range of 0.32 to 0.75 and a maximum occurred at aw of 0.75."

### Table V. "Regressions of ln reaction rates (ppm/hr) for the formation of pyrazine and 2-methylpyrazine at 95 C in NFDM as a function of water activity" (the three aw <= 0.75 points)

| compound | regression | n | R2 |
|---|---|---:|---:|
| pyrazine | ln y = 3.318 x - 2.222 | 3 | 0.986 |
| 2-methylpyrazine | ln y = 2.806 x - 2.174 | 3 | 0.990 |

Check: exp(3.318 x 0.319 - 2.222) = 0.312 vs Table IV 0.324; at 0.747: 1.29 vs 1.371. Consistent
with a fit to the measured aw values.

## 4. Kinetic numbers the repository can use

Registry mapping: 2-methylpyrazine -> `methylpyrazine`; 2,5-dimethylpyrazine -> `2_5_dimethylpyrazine`;
2,3-dimethylpyrazine -> `2_3_dimethylpyrazine`; **pyrazine (parent) -> not in registry as a molecule**
(class row `pyrazines` only); 2-methoxypyrazine (IS) -> class `methoxypyrazines`; lysine, glucose ->
not in `compounds.yml` (`reactive_lysine` is a marker).

Shared conditions for the solution rows: 100 mmol/L L-lysine·HCl + 100 mmol/L D-glucose, 0.1 M
buffer (citrate-phosphate at pH 5.0 and 7.0; borate at pH 9.0), initial pH as stated (drift
unquantified), capped Pyrex tubes in a water bath, absolute IS quantification (purge-and-trap NPD),
k = slope of concentration vs time (pseudo-zero order in product), 12-22 points over up to 24 h.

| quantity | value | unit | conditions | reaction order | source location | evidence class |
|---|---|---|---|---|---|---|
| pyrazine rate, pH 7.0 | 1.346 / 0.159 / 0.0957 (= 0.280 / 0.0331 / 0.0199 µmol L-1 min-1) | ppm/h | 95 / 85 / 75 C | pseudo-zero | Table I, p. 200 | measured_rate |
| pyrazine rate, pH 5.0 | 0.0938 / 0.0232 / 0.00356 (0.0195 / 0.00483 / 0.000741) | ppm/h | same | pseudo-zero | Table I | measured_rate |
| pyrazine rate, pH 9.0 | 3.596 / 0.490 / 0.214 | ppm/h | same | pseudo-zero | Table I (= ch. 7) | measured_rate |
| methylpyrazine rate, pH 7.0 | 1.367 / 0.276 / 0.0945 (0.242 / 0.0489 / 0.0167) | ppm/h | same | pseudo-zero | Table I | measured_rate |
| methylpyrazine rate, pH 5.0 | 0.00636 / 0.00279 (95 / 85 C only) | ppm/h | r2 0.890 / 0.912 | pseudo-zero | Table I | measured_rate (weak) |
| 2,5-dimethylpyrazine rate, pH 7.0 | 0.0630 / 0.00949 / 0.00209 (0.00971 / 0.00146 / 0.000322) | ppm/h | same | pseudo-zero | Table I | measured_rate |
| 2,5-DMP and 2,3-DMP at pH 5.0; 2,3-DMP at pH 7.0 | below quantification | — | 75-95 C | — | text p. 202 | level_only (null) |
| Ea, pyrazine, pH 7.0 / 5.0 / 9.0 | 33.4 / 41.8 / 35.8 (139.7 / 174.9 / 149.8) | kcal/mol (kJ/mol) | 75-95 C, 3 points | Arrhenius on zero-order k | Table II, p. 201 | measured_barrier |
| Ea, methylpyrazine, pH 7.0 / 9.0 | 34.0 / 36.7 (142.3 / 153.6) | kcal/mol (kJ/mol) | same | same | Table II | measured_barrier |
| Ea, 2,5-dimethylpyrazine, pH 7.0 / 9.0 | 43.7 / 41.9 (182.8 / 175.3) | kcal/mol (kJ/mol) | same | same | Table II | measured_barrier |
| Ea, 2,3-dimethylpyrazine, pH 9.0 | 44.8 (187.4) | kcal/mol (kJ/mol) | same | same | Table II | measured_barrier |
| k(pH 9)/k(pH 7), pyrazine / methylpyrazine / 2,5-DMP | 2.2-3.1 / 1.5-2.1 / 2.6-3.2 across 75-95 C | — | same system | — | derived from Table I | within-study ratio |
| k(pH 9)/k(pH 5), pyrazine / methylpyrazine | 21-60 / 151-446 | — | same | — | derived from Table I | within-study ratio (buffer changes between pH 7 and 9) |
| total pyrazines, 2 h at 95 C, pH 5 / 7 / 9 | 0.0239 / 6.19 / 13.1 | ppm | same | — | Table III, p. 202 | level_only (end-of-cook; validation) |
| distribution, 2 h at 95 C | Table III percentages | % | same | — | Table III | within-study ratio |
| NFDM pyrazine rate at aw 0.319 / 0.583 / 0.747 / 0.841, 95 C | 0.324 / 0.681 / 1.371 / 1.318 | "ppm/h" in standardized units per initial dry weight (not a concentration) | 5 g powder, 0.5-3 h | pseudo-zero | Table IV, p. 203 | measured_rate, unit non-transferable; use only as ratios (1 : 2.10 : 4.23 : 4.07) |
| NFDM methylpyrazine rate, same aw | 0.271 / 0.626 / 0.886 / 0.718 | same | same | pseudo-zero | Table IV | as above (1 : 2.31 : 3.27 : 2.65) |
| ln k vs aw slope, 0.32-0.75, pyrazine / methylpyrazine | 3.318 / 2.806 per unit aw | — | 95 C NFDM | — | Table V, p. 206 | within-study shape |
| pH-vs-k regression R2 0.974 / 0.999; Figures 1-4 | — | — | — | — | Figures | figure_only |

## 5. Flags

1. **Buffer confound in the pH ladder.** pH 5 and 7 are citrate-phosphate, pH 9 is borate; phosphate
   and citrate are Maillard catalysts, borate complexes sugars. The pH 7 -> 9 step therefore changes
   two things. The authors' assertion that buffer type made no difference is unsupported by printed
   data.
2. **Table IV/V units are not concentrations.** External-standard response factors applied to a
   directly purged dry powder, then normalised to initial dry weight and multiplied by 100. The
   numbers are a purge-recovery-weighted index whose absolute scale cannot be reconstructed; only
   the ratios across aw carry information, and even those assume purge efficiency from the powder is
   aw-independent (it is not obviously so: moisture affects desorption from a powder).
3. **Re-print discrepancies vs ch. 7** (same lysine-glucose pH 9.0 data): 2,5-DMP 95 C r2 0.976 here
   vs 0.995 there; 2,3-DMP intercepts -0.0057 / -0.00017 / -0.00086 here vs -0.00569 / -0.00173 /
   -0.00860 there; Ea 2-methylpyrazine 36.7 vs 36.6; Ea 2,5-DMP 41.9 vs 42.3. My refit of the
   (identical) k values reproduces ch. 7's Ea (36.6, 42.3), so ch. 7 is the version to carry for
   those two barriers.
4. **Three-point Arrhenius, curved.** Same caveat as ch. 7 (pairwise Ea for pH 9 pyrazine 218 vs 86
   kJ/mol). The "minimum at pH 7" in Ea (33.4 vs 35.8 and 41.8 kcal/mol) rests on three points per
   pH and no CIs; the authors themselves decline to interpret it further.
5. **Weak rows**: pyrazine pH 7 / 85 C (r2 0.899, intercept 0.609 ppm — a large intercept for a k
   of 0.159 ppm/h); 2,5-DMP pH 7 / 75 C (n 14, r2 0.857); methylpyrazine pH 5 (r2 0.89-0.91). The
   pH 7 / 85 C pyrazine point breaks the monotone k(pH 9)/k(pH 7) pattern (3.08 vs 2.67 and 2.24).
6. **pseudo-zero order = initial rate** (0.3 % conversion at the most productive condition); rates
   must be stored with [lysine] = [glucose] = 0.1 M attached; no order in reactants is determined.
7. **"ppm" volume basis** (10 mL reaction -> 15 mL analysed) unstated, as in ch. 7; affects absolute
   k by a possible factor 1.5, not ratios or Ea.
8. **pH drift during heating unquantified** at all three pH values; formic/acetic acid production
   would push the pH 7 and 9 systems down, more so at 95 C and 24 h.
9. Lysine used as the monohydrochloride: 0.1 M chloride present; at pH 9.0 the epsilon-amine (pKa
   ~10.5) is mostly protonated — the "two amino groups" rationale is only partly realised.
10. Reference 18 (Part One) is cited as ACS Advances in Chemistry "NO. 388 ... p. 76"; the chapter is
    ACS Symposium Series 388, pp. 76-91.
11. Not measured: reactant loss, dicarbonyls, browning, Strecker aldehydes; no temperature at or
    above 100 C; no alanine or glycine; no 2,6-dimethylpyrazine; NFDM sampling times (five, 0.5-3 h)
    not individually printed; NFDM aw after heating not reported.
