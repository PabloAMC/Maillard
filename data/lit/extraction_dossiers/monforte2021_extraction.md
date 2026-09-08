# Monforte 2021 — EXTRACTION (phenylalanine 2.4 mmol/L + gallic acid 2.4 mmol/L, ± glucose 2.4 mmol/L, ± Fe(II) 7.5 mg/L and Cu(II) 0.4 mg/L, 10 mmol/L tartrate pH 3.4, headspace swept with oxygen; phenylacetaldehyde release followed on-line by PTR-TOF-MS at 40-80 C; Weibull rate constants and a 73.5 kJ/mol Arrhenius barrier)
### A five-temperature barrier for phenylalanine -> phenylacetaldehyde, but for the wine (o-quinone / metal / oxygen) Strecker route at pH 3.4, measured as headspace release, with glucose as an inhibitor — not the sugar-dicarbonyl Strecker step R07 models.

**Source on disk:** `data/articles/Monforte2021.pdf` (6 pp., owner's download, 2026-09-08). Read from
the text layer (`scratchpad/articles/Monforte2021.txt`); Tables 1 and 2 came through clean and are
re-typed below. Figure 2B (bar chart of Weibull k for four substrate/metal combinations) is not read;
the two k values the text quotes from it are recorded as text values. Figures 3 and 4 (release
traces) not read. No Supporting Information is referenced.

## 0. Identity

| field | value |
|---|---|
| Title | "Phenylacetaldehyde real-time release kinetics in wine like model solutions" |
| Authors | Ana Rita Monforte, Sara I.F.S. Martins, António César Silva Ferreira* (Universidade Católica Portuguesa, Porto; Wageningen University; Stellenbosch; Cork Supply Portugal) |
| Venue | Food Chemistry 364 (2021) 128948. Received 21 July 2020, revised 10 December 2020, accepted 22 December 2020, online 1 January 2021 |
| DOI | 10.1016/j.foodchem.2020.128948 |
| Naming | PA = phenylacetaldehyde; SA = Strecker aldehyde(s); SD = Strecker degradation; "oxidation model" = phenylalanine + gallic acid + Fe2+ + Cu2+; "release" = ppbv in the swept headspace per unit time; "cumulative" = the time integral of release |
| Companions | Monforte, Martins & Silva Ferreira 2018 (JAFC 66, 2459; the GC-based precursor study, not on disk); Chan & Reineccius 2005 (Ea 21.5 kJ/mol for PA from the Maillard reaction at pH 7, cited not read); da Costa 2004 (76.7 kJ/mol in lager beer, cited not read); Hofmann, Münch & Schieberle 2000 (`hofmann2000_extraction.md`, the pH dependence cited here) |

## 1. Why it matters

Programme 6 (`tasks/roadmap_for_scientists.md` 5c) needs per-amino-acid Strecker rates or barriers at
two or more temperatures; phenylalanine -> phenylacetaldehyde is one of the six missing Strecker
rows. This paper prints a rate constant at five temperatures (40-80 C) with an Arrhenius barrier of
73.5 ± 0.8 kJ/mol. Three things limit its transport into the sugar-path lane: (i) the oxidant is the
o-quinone of gallic acid formed with Fe/Cu and a headspace swept with oxygen, i.e. the wine
oxidation route, not glyoxal or methylglyoxal from a sugar — and glucose *lowers* the rate here;
(ii) the measured quantity is the cumulative headspace release of PA (ppbv integrated over time),
fitted with a Weibull function, so k is a release rate lumping formation, partition and purge, not a
liquid-phase formation rate; (iii) pH 3.4 tartrate. What it does give the repository: a measured
barrier for the quinone-mediated Strecker of Phe (a candidate for a phenolic-rich matrix, or as a
bound on the acid-side behaviour), the pH ordering 7 > 5 > 3.4 at 80 C, a four-fold metal effect,
a six-fold SO2 suppression and, for the discussion of oxidant-limited channels (B17), an example of
the opposite regime: oxygen in excess.

## 2. Methods as they matter to a model

- **Model solutions.** "equimolar amounts of phenylalanine (2.4 mM, 0.4 g/L), gallic acid (2.4 mM,
  0.4 g/L) and glucose (2.4 mM, 0.4 g/L) were added to a tartrate buffer at pH 3.4. Tartrate buffer
  (10 mM, pH = 2.9) was prepared by mixing 7.5 mM of tartaric acid with 2.5 sodium tartrate
  dihydrate and the pH was adjusted to 3.4 with NaOH 1 M." Check: 0.4 g/L / 165.19 = 2.42 mmol/L
  Phe; 0.4 / 170.12 = 2.35 mmol/L gallic acid; 0.4 / 180.16 = 2.22 mmol/L glucose (the "2.4 mM" is
  nominal for the last two). Metals: **Cu2+ 0.4 mg/L = 6.3 µmol/L** (CuSO4·5H2O), **Fe2+ 7.5 mg/L =
  134 µmol/L** (FeSO4·7H2O), "as these concentrations are normally found in white wines". **SO2 50
  mg/L = 0.78 mmol/L** (form of addition not stated). pH 5 and 7 solutions: how the pH was set is not
  stated.
- **Vessel and flow.** 50 mL of solution in a 100-mL three-neck round-bottom flask in a silicone
  bath; "the headspace (50 mL) was swept with oxygen"; neck A gas in, neck C to the PTR-MS inlet,
  neck B for additions. **The oxygen flow rate is not stated.** Runs of 60 min.
- **Detector.** PTR-TOF-MS 8000 (Ionicon), inlet 125 C, drift tube 90 C, 3.8 mbar, 900 V, E/N 131
  Td; PA followed at m/z 121.065 (C8H9O+); concentrations in ppbv from primary-ion counts with a
  common rate coefficient k = 2e-9 cm3/s (Deuscher 2019) — i.e. **no compound-specific
  calibration**; linearity checked by injecting "measured amounts" of PA (values not given).
  Three replicates per condition, spectra averaged. No liquid-phase concentration of PA is ever
  reported.
- **Kinetic treatment.** "Release kinetics was analyzed by converting the curves of release to
  cumulative release vs. time. The cumulative release profiles were fitted to the Weibull model
  (Mateus, Lindinger, Gumy, & Liardon, 2007), using Excel for the calculation of the process rate
  constant k (min-1) and the shape parameters (n)." The Weibull equation itself is **not printed**;
  Mateus 2007 uses C(t) = C_inf [1 - exp(-(k t)^n)] (from memory of that paper, unverified here). n
  "indicates concavity or convexity ... when it takes values below and above 1"; Solver
  least-squares; R2 reported. Arrhenius: "activation energy of 73.5 ± 0.8 kJ/mol was obtained with a
  coefficient of determination of 0.999" on the five k of Table 1.
- **Which pot carries the temperature series.** Section 2.2: "Temperature effects were studied in the
  oxidation model (phenylalanine + gallic acid + Fe2+ & Cu2+) at 40, 50, 60, 70 and 80 C. Effect of
  pH was determined in the same model at 80 C for pH 3.4, 5 and 7." Section 3.2: "solution with
  phenylalanine, glucose, gallic acid and metals was further studied at three pH ... and five
  temperatures". **The two statements disagree on whether glucose was present** in the Table 1
  runs (see Flags 3).
- **Units.** k in min-1 as printed; 1 min-1 = 1.667e-2 s-1. Release values (Table 2, Fig. 3 numbers)
  are "ppb" of PA in the swept headspace; they cannot be converted to mol/L without the flow rate
  and the partition coefficient, neither given.

## 3. Tables re-typed

### Table 1. "Rate of formation (k) and shape parameter (n) from the fitting of the Weibull model for the three pH (3.4, 5 and 7) and the five temperatures (40, 50, 60, 70 and 80 C)"

Pot: phenylalanine + gallic acid + metals (+ glucose per section 3.2; see Flags 3), 60-min runs,
oxygen-swept headspace.

| series | condition | k (min-1) | n | R2 |
|---|---|---:|---:|---:|
| pH (T = 80 C) | pH 3.4 | 0.30 ± 0.02 | 0.56 ± 0.04 | 0.997 |
| | pH 5 | 0.32 ± 0.01 | 0.56 ± 0.03 | 0.999 |
| | pH 7 | 0.36 ± 0.05 | 0.56 ± 0.03 | 0.998 |
| Temperature (pH 3.4) | 40 C | 0.014 ± 0.008 | 0.11 ± 0.03 | 0.999 |
| | 50 C | 0.033 ± 0.004 | 0.24 ± 0.06 | 0.999 |
| | 60 C | 0.080 ± 0.003 | 0.56 ± 0.04 | 0.999 |
| | 70 C | 0.174 ± 0.01 | 0.56 ± 0.05 | 0.999 |
| | 80 C | 0.354 ± 0.02 | 0.56 ± 0.04 | 0.999 |

Note the 80 C / pH 3.4 entry appears twice (0.30 ± 0.02 in the pH series, 0.354 ± 0.02 in the
temperature series): two runs of nominally the same condition, 18 % apart.

**Arrhenius arithmetic (mine).** Unweighted least squares on the five (1/T, ln k) pairs: Ea = 74.7
kJ/mol, ln A = 24.43 (A = 4.1e10 min-1), R2 = 0.9997; two-point 40-80 C: 74.2 kJ/mol. The printed
73.5 ± 0.8 (R2 0.999) is reproduced within 1.2 kJ/mol; carry the printed value. k(80)/k(40) = 25.3
("increases 25 times").

### Table 2. "Phenylacetaldehyde release values (ppb) before and after SO2"

Pot: phenylalanine + gallic acid + glucose + metals; SO2 50 mg/L added at 20 min (the text says 25
min in one place); release = ppbv per unit time (the unit of time is not stated).

| series | condition | before SO2 (t = 20 min) | after SO2 (t = 25 min) | end (t = 60 min) |
|---|---|---:|---:|---:|
| pH (T = 80 C) | pH 3.4 | 1.4E-02 | 9.4E-04 | 2.2E-03 |
| | pH 5 | 1.6E-02 | 1.1E-03 | 2.4E-03 |
| | pH 7 | 2.3E-02 | 1.8E-03 | 2.7E-03 |
| Temperature (pH 3.4) | 60 C | 6.5E-04 | 1.5E-04 | 1.2E-04 |
| | 70 C | 2.8E-03 | 7.3E-04 | 4.6E-04 |
| | 80 C | 1.4E-02 | 9.9E-04 | 2.7E-03 |

Ratios before/after (mine): pH 3.4 14.9, pH 5 14.5, pH 7 12.8 ("15 times ... 13 times"); 60 C 4.3,
70 C 3.8, 80 C 14.1 ("4 times for 60 and 70 C and 14 times for 80 C"). End/after: 80 C 2.7 (release
recovers), 60 C 0.8, 70 C 0.6 ("for 60 and 70 C the release is constant while for 80 C an increase").

### Numbers in the running text

- Fig. 2B values quoted: Phe + gallic acid, 80 C, pH 3.4: **k = 0.1 min-1 without metals, 0.42
  min-1 with Fe2+ + Cu2+** ("increases four times with metals addition (0.1 to 0.42)"). Glucose
  addition "decreases the rate of PA formation" (value not quoted; Fig. 2B, FIGURE-ONLY). Weibull
  fits r2 > 0.95 for all Fig. 2 models.
- Sequential additions (Fig. 3; units not stated, presumably the same "release" as Table 2): Fe
  added after glucose 6.66e-5 ± 2.24e-6 vs Fe added before glucose 8.85e-6 ± 2.48e-7 (p < 0.05);
  Cu: 3.14e-5 ± 4.95e-7 (glucose before) vs 1.14e-4 ± 7.07e-7 (glucose after); end of reaction (60
  min) with metals after glucose 1.08e-4 ± 4.24e-7.
- SO2 from t = 0 (pH 3.4, 80 C, full pot): release 0.0375 without vs 0.0058 with SO2 ("6 times
  lower"; units not stated).
- Shape parameter: "approximately one, meaning that the curve fitting is similar to the first-order
  kinetics" (text) — Table 1 prints 0.11-0.56 (Flags 2).

## 4. Kinetic numbers the repository can use

Registry mapping: phenylacetaldehyde -> `phenylacetaldehyde`; phenylalanine, gallic acid, glucose,
sulfur dioxide / bisulfite, Fe(II), Cu(II), the PA-bisulfite adduct -> not in registry.

| quantity | value | unit | conditions | reaction order | source location | evidence class |
|---|---|---|---|---|---|---|
| Weibull release rate constant k for PA, 40 / 50 / 60 / 70 / 80 C | 0.014 / 0.033 / 0.080 / 0.174 / 0.354 (SD 0.008 / 0.004 / 0.003 / 0.01 / 0.02) | min-1 | Phe 2.4 mmol/L + gallic acid 2.4 mmol/L + Fe(II) 134 µmol/L + Cu(II) 6.3 µmol/L (+ glucose 2.4 mmol/L per 3.2), 10 mmol/L tartrate pH 3.4, 50 mL, O2-swept headspace, 60 min; n = 0.11 / 0.24 / 0.56 / 0.56 / 0.56 | Weibull on cumulative headspace release (n as printed; not a chemical order) | Table 1, p. 4 | measured_rate (of release, lumped; see Flags 1-2) |
| Ea for k | 73.5 ± 0.8 (R2 0.999); my refit 74.7 | kJ/mol | same, 40-80 C, five points | Arrhenius on the Weibull k | text 3.2, p. 4 | measured_barrier (quinone/metal/O2 Strecker route at pH 3.4, not the sugar-dicarbonyl route) |
| k at pH 3.4 / 5 / 7, 80 C | 0.30 / 0.32 / 0.36 (n = 0.56 all) | min-1 | same pot | Weibull | Table 1 | measured_rate; as ratios 1 : 1.07 : 1.20 within_study_ratio |
| metal effect on k, Phe + gallic acid, 80 C, pH 3.4 | 0.1 -> 0.42 (x 4.2) | min-1 | no glucose | Weibull | text 3.1 quoting Fig. 2B | within_study_ratio (the two values are printed in text; the bars are figure-only) |
| glucose effect on k | "decreases" | — | 80 C | — | Fig. 2B | figure_only |
| SO2 (50 mg/L) suppression of release | 0.0375 -> 0.0058 (x 6.5 lower) | release, unit not stated | pH 3.4, 80 C, SO2 from t = 0 | — | text 3.3 | within_study_ratio |
| release drop on SO2 addition at 20-25 min | x 14.9 / 14.5 / 12.8 (pH 3.4 / 5 / 7, 80 C); x 4.3 / 3.8 / 14.1 (60 / 70 / 80 C, pH 3.4) | — | same pot | — | Table 2 | within_study_ratio |
| release values before / after SO2 / end | Table 2 | ppb (headspace) | as above | — | Table 2 | level_only (headspace, no flow rate; not convertible to mol/L) |
| sequential-addition release values | 6.66e-5, 8.85e-6, 3.14e-5, 1.14e-4, 1.08e-4 | not stated | Fig. 3 protocol | — | text 3.1 | level_only (units missing) |
| cumulative PA profiles, all pots | — | ppb | — | — | Figs. 2A, 3, 4 | figure_only |

**Barrier context inside the corpus.** The roadmap carries Cremer & Eichner 2000 (Ea 115-124 kJ/mol
for Strecker aldehydes, cited not read) and the paper itself cites 21.5 kJ/mol (Chan & Reineccius
2005, Maillard, pH 7) and 76.7 kJ/mol (beer). The 73.5 kJ/mol here is for a different oxidant
(o-quinone) at a different pH; it should be stored under the phenolic/quinone route, not as the
R07 barrier for phenylalanine on glyoxal or methylglyoxal.

## 5. Flags

1. **Release, not formation; headspace, not liquid.** PA is measured as ppbv in an oxygen stream
   sweeping a 50-mL headspace at an unstated flow; k is fitted on the time integral of that signal.
   It lumps chemical formation, liquid-gas partition (temperature-dependent on its own: PA's Henry
   constant changes about threefold between 40 and 80 C) and purge dynamics. Part of the 73.5
   kJ/mol may be partition, not chemistry. No liquid concentration of PA is ever given, so no yield
   per mole of phenylalanine can be computed.
2. **The Weibull shape parameter is 0.11-0.56, not 1.** The text's "approximately one, meaning ...
   first-order" contradicts Table 1. With n = 0.11 (40 C) and 0.24 (50 C) the low-temperature k are
   not on the same footing as the n = 0.56 points, and the Arrhenius line mixes them; Weibull k
   values with different n are not directly comparable rate constants. The equation is not printed.
3. **Glucose in or out of the Table 1 pot?** Section 2.2 says the temperature and pH series were in
   the oxidation model (no glucose); section 3.2 says the pot had glucose. Since glucose lowers the
   rate (Fig. 2B), the identity of the pot matters for any reuse; unresolved from the paper.
4. **Oxidant in excess.** Oxygen-swept headspace, Fe(II) 134 µmol/L, Cu(II) 6.3 µmol/L, gallic acid
   2.4 mmol/L: the Strecker oxidant here is the o-quinone continuously regenerated by O2 and metals.
   This is the opposite of the oxidant-limited pots B17 diagnosed (`kinetic_core_b17_prereg.md`
   section 6); the rates are not transferable to a sealed sugar pot.
5. **pH 5 and 7 buffers not described**; tartrate at 10 mmol/L has little capacity at pH 7.
6. **PTR-MS quantification** uses a generic rate coefficient (2e-9 cm3/s) and no PA-specific
   calibration factor; ppbv values carry a systematic scale uncertainty (typically ±30 %). Ratios
   within the study are unaffected.
7. **Duplicate 80 C / pH 3.4 entries** (0.30 vs 0.354 min-1) show the run-to-run spread (18 %) is
   larger than the printed SDs (0.02).
8. **60-min runs**; at 40 C with k = 0.014 min-1 and n = 0.11 the cumulative curve is nowhere near
   complete, so the asymptote and k are extrapolations.
9. **Units missing** for the Fig. 3 release numbers and for the SO2 comparison (0.0375 vs 0.0058).
10. **The bisulfite-PA adduct** (hydroxysulfonate) reversibility at 80 C (release recovers from
    9.9e-4 to 2.7e-3 between 25 and 60 min) is a matrix effect for wine/beer, not for the meat
    analogue matrices the repository models; recorded for completeness.
11. **Registry**: only phenylacetaldehyde is keyed; phenylalanine is not a species in
    `compounds.yml`.
