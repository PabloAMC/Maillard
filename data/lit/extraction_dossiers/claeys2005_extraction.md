# Claeys, De Vleeschouwer & Hendrickx 2005 — EXTRACTION (equimolar L-asparagine + D-glucose 0.01 mol/L in 0.05 M citrate pH 6, plus 0.01 mol/L glutamine / cysteine / lysine / alanine, sealed inox tubes, 140-200 C; two consecutive first-order reactions fitted to the acrylamide time course; k_F, k_E, Ea_F, Ea_E ± SE at T_ref = 160 C for five systems)
### The dilute-aqueous half of the shipped lane's fit corpus, read from the primary table for the first time: five acrylamide FORMATION and five ELIMINATION constants at one reference temperature, an elimination barrier of 167.2 ± 4.3 kJ/mol that its own laboratory's later paper contradicts, and a competitor panel whose whole content is five (k_E/k_F) ratios.

**Source on disk:** `data/articles/claeys2005.pdf` (6 pp., Biotechnol. Prog. 21 (2005) 1525-1530).
Read from the `pdftotext -layout` text layer (`scratchpad/articles/claeys2005.txt`); Tables 1 and 2
came through clean and are re-typed in full below. Figures 1 (relative acrylamide vs time at
160 C), 2a-e (acrylamide vs time at four temperatures for the five systems, with the fitted curves)
and 3 (the Arrhenius plot of k_F and k_E for all five systems) are images: **every acrylamide
concentration in this paper is figure-only, and there is not one absolute acrylamide level printed
anywhere.** No supplementary material. Repo status before this dossier: Claeys 2005 is named in
`src/kinetic_core/parameters_acrylamide.py` (the MEASURED BACKBONE header, policy 3, GROUP 1 and
GROUP 2 comments, and the fitted-step conditions string) and is declared FIT in
`FIT_HOLDOUT_DECLARATION.md` D.5, but has **no extraction dossier** — the numbers reach the
repository through `k3_final_parameter_inventory.md` only.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of Amino Acids on Acrylamide Formation and Elimination Kinetics" |
| Authors | Wendie L. Claeys, Kristel De Vleeschouwer, Marc E. Hendrickx (corresponding) — Laboratory of Food Technology, Department of Food and Microbial Technology, Faculty of Applied Bioscience and Engineering, Katholieke Universiteit Leuven, Kasteelpark Arenberg 22, B-3001 Heverlee, Belgium |
| Venue | Biotechnol. Prog. 2005, 21 (5), 1525-1530. Accepted 8 August 2005, web 20 September 2005 |
| DOI / article ID | 10.1021/bp050194s (printed as `BP050194S`) |
| Naming | AA = acrylamide; R = "Reactants", a single lumped species with C_R0 = 0.01 M; D = the acrylamide sink, written "AA-protein complex / AA degradation product / ..."; k_F, k_E = first-order formation and elimination constants; Ea_F, Ea_E their barriers |
| Companion (cited, NOT on disk) | Claeys, De Vleeschouwer & Hendrickx (2005b), "Kinetics of acrylamide formation and elimination during heating of an asparagine-sugar model system", JAFC 53:9999-10005 — reference (19), the paper that introduced this two-consecutive-first-order scheme. **The repository does not have it.** |
| Companions on disk | `devleeschouwer2006_extraction.md` (pH 4/6/8), `devleeschouwer2007_extraction.md` (a_w 0.34-0.92), `devleeschouwer2008_extraction.md` (a_w 0.88-0.99), `devleeschouwer2009_extraction.md` (Part II, competitors, low moisture), and `knol2005_extraction.md` / `knol2010_extraction.md` (the Wageningen second laboratory) |

## 1. Why it matters

The shipped acrylamide lane is wave B3; its constants are in
`src/kinetic_core/parameters_acrylamide.py` and its network in `src/kinetic_core/acrylamide.py`.
**This paper is one of the three FIT sources that lane was built on**, and it is the *only* one of
the three in dilute aqueous solution — the registry's own policy 4 says so in as many words: "Claeys
is dilute aqueous (a_w ~ 1.0), De Vleeschouwer is a freeze-dried powder at a_w 0.92, and the
extrusion benchmark is a_w 0.35." The fitted steps' declared water activity,
`FITTED_AW_OF_DECLARATION = 1.0`, is this paper's water activity, because "the panel rows that
actually identify them" are these.

Concretely, this paper is where the repository gets four things:

1. **The acrylamide ELIMINATION constant of the Leuven laboratory in water**: k_Eref = 111.1e-3
   min^-1 at 160 C with Ea_E = 167.21 ± 4.30 kJ/mol. The lane's GROUP 2 comment quotes 0.1111
   Claeys / 0.100 De Vleeschouwer / 0.0881 Knol and 167.2 / 113.2 / 85.1 — **the first two of those
   six numbers are the control row of Table 2 here**, and both are confirmed against the primary
   table below. The lane's own note records that the retired FAST lane had no elimination step at
   all and was ~40x under-responsive to time; the panel's remaining complaint is that the time
   *shape* is inverted against a measurement (28 -> 912 -> 1459 ppb over 10-30 min in a 0.5 %
   glucose + 0.5 % asparagine solution). Section 4 does the arithmetic on what this paper's
   elimination constant implies for a 30-minute window.
2. **The Int1 partition** — indirectly. `k_int1_mel` is fitted, not carried, and the registry says
   why: De Vleeschouwer's k_M is marked "NO PHYSICAL MEANING" by its own authors, so the partition
   is "fitted here against Claeys' dilute-aqueous lumped constants, which is the only place in the
   corpus where the partition is observable." That claim rests on this paper's k_F, and section 4
   says exactly what k_F is and is not.
3. **The competitor panel.** The four fitted competitor channels (`k_gln_glc`, `k_lys_glc`,
   `k_ala_glc` and their acrylamide-scavenging partners) are fitted against this paper's four
   competitor systems. The registry's expectation "Claeys' alanine row is statistically
   indistinguishable from the control, so this constant is expected to come out at the bottom of
   its range" is confirmed below from Table 2's significance letters.
4. **The refusal of a pH term.** Policy 3 in the registry states "Claeys is 0.05 M citrate at
   pH 6"; confirmed. This paper varies amino acid, not pH.

What it is **not**: it is not a second laboratory. Claeys, De Vleeschouwer and Hendrickx are the
same Leuven group as De Vleeschouwer 2006 / 2007 / 2008 / 2009 I / 2009 II. Reading this dossier
next to `knol2005_extraction.md` and `knol2010_extraction.md` is the point: those two are the
outside opinion, this one is the incumbent.

## 2. Methods as they matter to a model

- **Pot.** "a model system consisting of 0.0 1M L-asparagine (>=99.5 %) and 0.01 M D-glucose
  (>=99 %) dissolved in a 0.05 M citrate buffer pH 6, to which 0.01 M L-glutamine (>=99.5 %),
  L-cysteine (>=99.5 %), L-lysine (>=98 %), or L-alanine (>=99.5 %) was added (Sigma-Aldrich, USA)."
  The typography "0.0 1M" is a printing artefact for 0.01 M. So: **asparagine 10 mmol/L, glucose
  10 mmol/L, second amino acid 10 mmol/L where present, citrate 50 mmol/L, pH 6.** This is
  **20x more dilute than Knol 2005** (0.2 mol/L each) and roughly **270x more dilute in molar terms
  than De Vleeschouwer's a_w 0.92 powder** (see `devleeschouwer2009_extraction.md` section 2).
  **Water activity ~1.0.** No pH re-measurement is reported; 50 mM citrate against 10 mM reactants
  is a genuinely buffered system, unlike Knol's.
- **Vessel and heating.** "Samples were heated in airtight reactor tubes (custom-made, inox,
  8 mm x 100 mm) to avoid as much as possible side phenomena, such as water evaporation and
  absorption of oil, which affect AA formation." Thermostated oil bath, **140-200 C** (Table 2 says
  "between 140 and 200 C"; Table 1 lists 140, 160, 180, 200), immediate ice-water cooling.
  **The temperature was logged inside the tube every 2 s with type-T thermocouples (Ellab TM 9616)
  and the kinetic fit is performed on the Euler integral of the recorded temperature-time
  profile** — i.e. the heat-up and cool-down are corrected for, unlike Knol 2005. This is the same
  inox tube (8 x 100 mm) and the same instrumentation used by the whole De Vleeschouwer series.
- **Times.** Not tabulated. Table 1's footnote says the relative yields are the "Average percentage
  of nine samples taken after different treatment times at the temperature in concern", so **nine
  time points per temperature**; Table 2 says "number of data points = 35" for each system, i.e.
  roughly nine points at four temperatures with one lost. The actual times appear only on the
  Figure 2 x-axes.
- **Acrylamide.** GC-MS, essentially Biedermann et al. 2002. Extraction without derivatisation,
  clean-up, 1 uL cool-on-column onto HP-INNOWax 30 m x 250 um x 0.25 um with a 0.5 m x 530 um
  deactivated fused-silica precolumn; quadrupole MS in **positive chemical ionisation**, He carrier,
  CH4 reagent gas (Agilent 5973 inert). **Acrylamide quantified at m/z 72 against methacrylamide
  (m/z 86) added at the start of sample preparation**, with **butyramide (m/z 88) added before
  injection as a second internal standard so that losses during analysis could be accounted for.**
  Two internal standards, one of them tracking extraction. No LOD/LOQ printed.
- **Kinetic scheme.** Two consecutive first-order reactions, taken from Claeys 2005b (ref. 19):
  `Reactants --k_F--> AA --k_E--> D` where D = "AA-protein complex / AA degradation product / ...".
  Rate laws as printed: dC_R/dt = -k_F C_R (eq 1); dC_AA/dt = k_F C_R - k_E C_AA (eq 2);
  dC_D/dt = k_E C_AA (eq 3). At t = 0, C_AA = C_D = 0 and **C_R0 = 0.01 M**.
  **The "reactants" are ONE lumped species at 0.01 M.** Asparagine and glucose are not separate
  variables and neither is measured; the second amino acid is not in the equations at all.
- **Temperature dependence.** k = k_ref exp[(Ea/R)(1/T_ref - 1/T)] (eq 4, as printed) with
  **T_ref = 160 C = 433.15 K** — the same reference the registry uses. Non-linear regression
  (Gauss-Newton, SAS v8) on the integral over the recorded temperature-time profile, fitted on the
  total data set. **Uncertainties are standard errors (SE), not HPD intervals** (contrast Knol and
  De Vleeschouwer). Significance letters in Table 2 are "based on 95 % asymptotic confidence
  intervals".
- **Goodness of fit.** pseudo-R^2 = 1 - SS_residual/SS_corrected (eq 5, Schabenberger 1998), plus
  observed-vs-predicted bias plots, residual plots, and the Shapiro-Wilk W for normality of
  residuals (Pr < W, rejected below 0.05).
- **Parameter correlation, stated by the authors.** "For all models studied, a high correlation was
  observed between k_Fref and k_Eref (0.94-0.98) and between Ea_F and Ea_E (0.74-0.92). This,
  however, does not give any information about the physical relationships between the parameters
  but only about the fitting process." **This is the single most important methodological caveat in
  the paper** and is carried into every row of section 4.

## 3. Tables re-typed

### Table 1. "Effect of Amino Acids on the AA Content (%) Formed in an Equimolar Asparagine-Glucose Model System (0.01 M, pH 6) Heated at Different Temperatures"

Values are percentages relative to the control at the same temperature. Footnote a: "Average
percentage of nine samples taken after different treatment times at the temperature in concern."

| T (C) | control | glutamine | cysteine | lysine | alanine |
|---:|---:|---:|---:|---:|---:|
| 140 | 100 (a) | 154.91 | 54.05 | 76.15 | 100.12 |
| 160 | 100 | 139.56 | 57.03 | 63.58 | 104.45 |
| 180 | 100 | 277.35 | 64.61 | 83.40 | 112.94 |
| 200 | 100 | 321.68 | 69.37 | 91.38 | 112.24 |

**No absolute acrylamide concentration is printed anywhere in this paper.** Table 1 is entirely
relative; the absolute levels behind it are on the Figure 2 y-axes only.

### Table 2. "Effect of Amino Acids on Kinetic Parameters Describing AA Formation/Elimination in an Equimolar Asparagine-Glucose Model System (0.01 M, pH 6) Heated between 140 and 200 C (Tref = 160 C)"

Units as printed: **k_Fref (x 10^-3 min^-1)** and **k_Eref (x 10^-3 min^-1)** — both first order.
Uncertainties are **standard errors**. Superscript letters a-c: "Values of the same parameter with a
different letter are significantly different based on 95 % asymptotic confidence intervals; number
of data points = 35; pseudo-R^2 = 1 - SS(residual)/SS(corrected total)."

| system | k_Fref (x10^-3 min^-1) | k_Eref (x10^-3 min^-1) | Ea_F (kJ/mol) | Ea_E (kJ/mol) | pseudo-R^2 | Pr < W |
|---|---|---|---|---|---:|---:|
| control | 0.451 ± 0.023 **a** | 111.1 ± 8.9 **a** | 168.25 ± 3.80 **a** | 167.21 ± 4.30 **a** | 0.975 | 0.934 |
| glutamine | 1.640 ± 0.416 **b** | 274.1 ± 81.6 **a,b** | 166.8 ± 14.4 **a,b** | 103.9 ± 17.5 **b** | 0.912 | 0.620 |
| cysteine | 0.501 ± 0.116 **a** | 268.7 ± 81.6 **a,b** | 206.3 ± 13.5 **b** | 180.0 ± 17.0 **a,c** | 0.765 | 0.599 |
| lysine | 0.587 ± 0.074 **a** | 280.2 ± 43.4 **b** | 179.3 ± 7.7 **a,b** | 140.0 ± 9.0 **b,c** | 0.950 | 0.572 |
| alanine | 0.465 ± 0.034 **a** | 103.1 ± 12.0 **a** | 173.3 ± 5.4 **a,b** | 169.7 ± 6.1 **a,c** | 0.956 | 0.096 |

Reading the letters: **k_Fref** — only glutamine differs from the control. **k_Eref** — only lysine
differs from the control. **Ea_F** — only cysteine differs from the control. **Ea_E** — glutamine
and lysine differ from the control; cysteine and alanine do not.

### Numbers printed in the running text

| quantity | value | where |
|---|---|---|
| k_Eref/k_Fref ordering | "increased from the systems to which glutamine or alanine was added, over the control model system, to the model systems that contained lysine or cysteine" | Abstract, Results |
| cysteine's k_Eref/k_Fref vs control | "twice as high as for the control system" | Results |
| glucose + glutamine alone (no asparagine) | "only negligibly small amounts of AA were measured ... heated between 140 and 200 C for 20 min" | Results (preliminary experiment, no number) |
| reported amino-acid stock impurities (Becalski et al., cited) | 0.3 % asparagine in an aspartic acid standard; 0.8 % cysteine in a glutamine standard | Results |
| SH vs NH2 reactivity toward conjugated vinyls (Friedman, cited) | 100-300x | Results |
| glutamine's Arrhenius crossing | the k_E curve of the glutamine system crosses the control's "around 180 C"; the k_F curves "run parallel" | Results, Fig. 3 |
| Rydberg et al.'s contrary result (cited) | 35 mM added amino acids to potato at 180 C / 25 min: glutamine and glycine -76 to -70 %, lysine -57 %, alanine -14 %; endogenous free asparagine 17 mM; pH not stated | Results |

### Arithmetic on the printed constants (all mine)

**1. k_Eref/k_Fref, the ratio the paper's own argument turns on.**

| system | k_Eref/k_Fref |
|---|---:|
| glutamine | 274.1 / 1.640 = **167** |
| alanine | 103.1 / 0.465 = **222** |
| control | 111.1 / 0.451 = **246** |
| lysine | 280.2 / 0.587 = **477** |
| cysteine | 268.7 / 0.501 = **536** |

This reproduces the abstract's ordering exactly (glutamine or alanine < control < lysine or
cysteine) and the "twice as high as for the control" claim for cysteine (536/246 = 2.18).

**2. Relative to control, at T_ref = 160 C.**

| system | k_F / k_F(control) | k_E / k_E(control) |
|---|---:|---:|
| glutamine | **3.64** | **2.47** |
| cysteine | 1.11 | **2.42** |
| lysine | 1.30 | **2.52** |
| alanine | 1.03 | 0.93 |

Alanine is within 7 % of the control on both constants — the registry's expectation that
`k_ala_glc` should "come out at the bottom of its range" is what this row says.

**3. The competitor effect does not reproduce Table 1 from Table 2 alone.** At 160 C Table 1 says
cysteine cuts the net acrylamide to 57 % and lysine to 64 %, yet Table 2 gives cysteine and lysine
k_F values statistically indistinguishable from the control and k_E values only ~2.5x higher.
In this scheme, with k_F << k_E (0.451e-3 against 111.1e-3, a factor 246), the peak acrylamide is
**(k_F/k_E) x C_R0** to within a percent, so a 2.5x rise in k_E at unchanged k_F cuts the peak
2.5x (mine). Table 1's time-averaged cut at 160 C is only to 57 % (1.75x). Same order, so the
picture is consistent, but the fit's k_F/k_E correlation of 0.94-0.98 means the split between
"less formed" and "more removed" is **not identified by these data**. The authors say the same
thing in different words: "The different actions of cysteine and lysine are, however, not
distinguishable in the kinetic model."

**4. What k_F hides.** k_F is first order in a single lumped "Reactants" species at 0.01 M. It is
therefore a **pseudo-constant that hides the glucose (and the asparagine)**. If one insists on
reading the initiation as bimolecular at these concentrations, k_F/[Glc] = 0.451e-3/0.01 =
**0.0451 M^-1 min^-1 (mine, an assumption the paper does not license)** — which is 15x below
De Vleeschouwer 2009 Part I's k_INTg (1.70 M^-1 min^-1) and 15x below Knol 2005's k1 (0.668
M^-1 min^-1, unit inferred). That gap is not a contradiction: k_F is the *whole* route from
reactants to acrylamide, including the partition that sends ~95 % of the condensed material to
browning, whereas k_INTg and k1 are the condensation alone. **The comparison that is licensed is
k_F against a lumped route, not against a condensation constant.**

**5. Acrylamide half-life from k_E (mine, first order, no re-formation).**
t_1/2 = ln2 / k_E at T_ref = 160 C = ln2/0.1111 = **6.24 min** (control). Rolling the printed
Arrhenius out, k_E(T) = 0.1111 exp[(167210/8.314)(1/433.15 - 1/T)]:

| T | k_E (min^-1) | t_1/2 (min) |
|---:|---:|---:|
| 120 C (**extrapolated below the window**) | 9.80e-4 | **708** |
| 140 C | 1.16e-2 | **59.8** |
| 160 C (= k_Eref) | 0.1111 | **6.24** |
| 180 C | 0.855 | **0.81** |
| 200 C | 5.58 | **0.124** |

The 167 kJ/mol barrier makes this an extremely steep sink: it is essentially inert over half an
hour at 120 C and gone in seconds at 200 C. Compare Knol 2005's measured k6 series, whose
barrier of 85.1 kJ/mol gives 87 / 24.7 / 7.9 / 2.8 / 1.07 min at the same five temperatures — the
two laboratories agree at 160 C by construction and diverge by 8x at 120 C and by 9x at 200 C.
**That divergence, not the 160 C rate, is what a second-laboratory comparison is actually about.**

**6. Ea_F ~ Ea_E in the control, to within their SEs.** 168.25 ± 3.80 against 167.21 ± 4.30. With
a fitted Ea_F/Ea_E correlation of 0.74-0.92 this near-equality is at least partly a fitting
artefact, and the authors say so. But it has a consequence the repository should note: **in the
control system the ratio k_E/k_F is nearly temperature-independent**, so the acrylamide peak
height barely moves with temperature and only the peak *time* does. That is the opposite of what
`k_acr_cys` produces in the shipped lane, where Ea_E2 = 51.3 is half the elimination barrier and
the two channels deliberately cross over in temperature.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** `acrylamide` is keyed. **Asparagine, glucose,
glutamine, cysteine, lysine, alanine, and the lumped "Reactants" and "D" species are NOT in the
registry** (checked: 75 ids, no sugars, no free amino acids). Shared conditions for every row:
10 mmol/L asparagine + 10 mmol/L glucose (+ 10 mmol/L second amino acid), **0.05 M citrate,
pH 6, dilute aqueous a_w ~ 1.0**, sealed inox tube 8 x 100 mm, oil bath, **140-200 C**, nine
sampling times per temperature, 35 data points per fit, temperature-time profile integrated,
**T_ref = 160 C = 433.15 K** (identical to `T_REF_A_K`).

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| Reactants -> acrylamide (LUMPED whole route) | k_Fref, control | **0.451 ± 0.023** | 10^-3 min^-1 | as above, no second amino acid | **first order in a single lumped "Reactants" at 0.01 M** — hides both glucose and asparagine | Table 2 p. 1528 | measured_rate |
| " | Ea_F, control | **168.25 ± 3.80** | kJ/mol | 140-200 C | — | Table 2 | measured_barrier |
| **acrylamide -> D (ELIMINATION)** | **k_Eref, control** | **111.1 ± 8.9** | **10^-3 min^-1** (= 0.1111 min^-1) | as above | **first order in acrylamide** | Table 2 | measured_rate — **the dilute-aqueous elimination anchor of the shipped lane** |
| " | **Ea_E, control** | **167.21 ± 4.30** | kJ/mol | 140-200 C | — | Table 2 | measured_barrier |
| Reactants -> acrylamide, + glutamine | k_Fref / Ea_F | 1.640 ± 0.416 / 166.8 ± 14.4 | 10^-3 min^-1 / kJ/mol | + 10 mmol/L Gln | first order, lumped | Table 2 | measured_rate (k_F **significantly** above control) |
| acrylamide -> D, + glutamine | k_Eref / Ea_E | 274.1 ± 81.6 / 103.9 ± 17.5 | 10^-3 min^-1 / kJ/mol | + 10 mmol/L Gln | first order | Table 2 | measured_rate (k_E n.s. vs control; Ea_E **significantly** below) |
| Reactants -> acrylamide, + cysteine | k_Fref / Ea_F | 0.501 ± 0.116 / 206.3 ± 13.5 | " | + 10 mmol/L Cys | first order, lumped | Table 2 | measured_rate (k_F n.s.; Ea_F **significantly** above control) |
| acrylamide -> D, + cysteine | k_Eref / Ea_E | 268.7 ± 81.6 / 180.0 ± 17.0 | " | + 10 mmol/L Cys | first order — **NOT resolved into a cysteine-dependent channel** | Table 2 | measured_rate (both n.s. vs control; see Flags 5) |
| Reactants -> acrylamide, + lysine | k_Fref / Ea_F | 0.587 ± 0.074 / 179.3 ± 7.7 | " | + 10 mmol/L Lys | first order, lumped | Table 2 | measured_rate (both n.s.) |
| acrylamide -> D, + lysine | k_Eref / Ea_E | 280.2 ± 43.4 / 140.0 ± 9.0 | " | + 10 mmol/L Lys | first order | Table 2 | measured_rate (k_E **significantly** above control; Ea_E **significantly** below) |
| Reactants -> acrylamide, + alanine | k_Fref / Ea_F | 0.465 ± 0.034 / 173.3 ± 5.4 | " | + 10 mmol/L Ala | first order, lumped | Table 2 | measured_rate (both n.s. — the null arm) |
| acrylamide -> D, + alanine | k_Eref / Ea_E | 103.1 ± 12.0 / 169.7 ± 6.1 | " | + 10 mmol/L Ala | first order | Table 2 | measured_rate (both n.s. — the null arm) |
| competitor effect on the elimination/formation balance | k_Eref/k_Fref | Gln 167, Ala 222, control 246, Lys 477, Cys 536 | — | 160 C | — | derived from Table 2 (mine) | within_study_ratio |
| competitor effect on each constant | k_F ratio / k_E ratio vs control | Gln 3.64 / 2.47; Cys 1.11 / 2.42; Lys 1.30 / 2.52; Ala 1.03 / 0.93 | — | 160 C | — | derived (mine) | within_study_ratio |
| net acrylamide, relative to control | Gln 154.91 / 139.56 / 277.35 / 321.68; Cys 54.05 / 57.03 / 64.61 / 69.37; Lys 76.15 / 63.58 / 83.40 / 91.38; Ala 100.12 / 104.45 / 112.94 / 112.24 | % of control | 140 / 160 / 180 / 200 C, averaged over nine times | — | Table 1 p. 1527 | within_study_ratio (time-averaged; **not** a level) |
| acrylamide half-life from k_E, control | t_1/2 | 59.8 / **6.24** / 0.81 / 0.124 (and 708 at 120 C, extrapolated below the window) | min at 140/160/180/200 C | as above | — | derived from Table 2 (mine) | derived_assumption |
| acrylamide time courses, absolute levels, sampling times | — | — | — | — | — | Figs. 1, 2a-e, 3 | **figure_only** |
| glucose + glutamine without asparagine | "negligibly small amounts" of acrylamide, 140-200 C, 20 min | — | — | — | — | Results | level_only (no number printed) |

### Can these be put on the same basis as the trunk's constants? Step by step.

**T_ref matches exactly (160 C), and the unit on both constants is a plain first-order min^-1**, so
no transport arithmetic is needed to place these next to `k_int1_acr` and `k_acr_dp`. What differs
is what the constants *mean*.

**(a) The elimination is directly comparable, and is the incumbent the second laboratory must be
checked against.** k_Eref = 0.1111 min^-1 at 160 C, Ea_E = 167.21 ± 4.30. At the same T_ref:
De Vleeschouwer 2009 Part I glucose gives 0.10 ± 0.04 min^-1 with Ea_E = 113.2 ± 32.3, and Knol
2005's k6 gives 0.0881 ± 0.025 min^-1 with Ea = 85.1 ± 14. **The three rates agree to 1.26x.
The barriers do not: Claeys' 167.21 ± 4.30 and De Vleeschouwer's 113.2 ± 32.3 do not overlap
(162.9 vs 145.5) even though they come from the same laboratory**, and Knol's 85.1 ± 14 is further
away still. Claeys' interval is by far the tightest of the three (SE 4.30 on a two-parameter fit
to 35 points), and it is also the one whose fit reports a k_F/k_E correlation of 0.94-0.98 — so
its tightness should not be read as its being the best-determined. **A single (k_ref, Ea) pair
fitted to all three, as `parameters_acrylamide.py` GROUP 2 does, is the honest treatment, and this
dossier supplies the primary numbers for two of the three.**

**(b) The formation constant is NOT comparable to `k_int1_acr` as a like-for-like.** `k_int1_acr` is
first order in an unmeasured intermediate INT1, downstream of a bimolecular condensation. Claeys'
k_F is first order in a lumped "Reactants" that starts at 0.01 M and stands for *the whole route*.
The two have the same unit and the same T_ref but different denominators. What Claeys' k_F does
license is the registry's actual use of it: as the observable against which the **partition**
`k_int1_mel` is fitted, because k_F is the only constant in the corpus that measures how much of
the reactant pool arrives at acrylamide in a system where the condensation constant is
independently known. **It is a composite, and it should keep the label.**

**(c) The competitor rows are within-study ratios and nothing more.** Neither Table 1 nor Table 2
resolves a mechanism. The cysteine row's k_E is *not* a cysteine-dependent second-order constant —
it is the same first-order k_E refitted in a pot that happens to contain cysteine, and its SE
(±81.6 on 268.7, 30 %) makes it statistically indistinguishable from the control. The measured
second-order acrylamide + cysteine constant the lane carries (`k_acr_cys` = 49.36 M^-1 min^-1,
Ea 51.3 ± 1.5) comes from De Vleeschouwer 2009 Part II, not from here. **What this paper adds is
the ORDERING and the four ratios in section 3.1-3.2**, which is precisely what fitted channels
with bounds that allow zero can be scored against.

**(d) What cannot be transported.** No pH variation (one pH, 6, citrate-buffered). No water
activity or moisture variation. No absolute concentration of anything, so **no benchmark row**.
No time-resolved data. No 120 C point — the window starts at 140 C, so any use of these constants
at 120 C is an extrapolation, and the repository's declared temperature range for the Claeys rows
should say 140-200 C, not 120-200 C (Flags 2).

## 5. Flags

1. **k_F and k_E are correlated at 0.94-0.98 and the authors say the correlation "does not give any
   information about the physical relationships between the parameters but only about the fitting
   process".** Ea_F and Ea_E are correlated at 0.74-0.92. In a two-consecutive-first-order model
   fitted to the *product* curve only, this is structural: the rise is set by k_F and the fall by
   k_E, and the peak height by their ratio, so the pair is much better determined than either
   member. **Any adoption should carry the pair or the ratio, never k_E alone.** This also means
   the tight ±4.30 on Ea_E must not be read as high precision on the elimination barrier.
2. **The temperature window is 140-200 C, not 120-200 C.** Table 2's title says "Heated between 140
   and 200 C" and Table 1 lists only those four. `parameters_acrylamide.py` declares its Claeys-
   anchored fitted rows over `temperature_range_c = (120.0, 200.0)` and its conditions string says
   "Claeys 2005 dilute aqueous pH 6 (140-200 C)" — the conditions string is right and the range
   tuple is inclusive of Knol's and De Vleeschouwer's 120 C points, which is defensible for a
   *joint* fit but should not be read as Claeys covering 120 C.
3. **Not one absolute acrylamide concentration is printed.** Table 1 is percentages of control;
   Figures 1-3 carry the only levels. This paper can supply constants but **cannot supply a
   benchmark row**, and the yield per mol asparagine is not derivable from it.
4. **"0.0 1M" in the Materials section** is a printing artefact for 0.01 M; every table title and
   the model's C_R0 confirm 0.01 M. Recorded so nobody reads it as 0.1 M.
5. **The cysteine and glutamine rows have SEs of 30 % and 30 % on k_E and pseudo-R^2 of 0.765 and
   0.912.** Cysteine's fit is the worst in the paper (0.765) — expected, since cysteine is the one
   system where the assumed mechanism (a first-order elimination) is wrong: De Vleeschouwer 2009
   Part II later shows the cysteine effect is a **second-order, cysteine-concentration-dependent**
   elimination. **The cysteine row of Table 2 is therefore a mis-specified fit, and its k_E and
   Ea_E should not be adopted as constants**; its value is as the earlier, cruder measurement that
   Part II supersedes.
6. **Alanine's Pr < W = 0.096** is the closest any system comes to failing the Shapiro-Wilk
   normality test at 0.05. Not a rejection, but the alanine residuals are the least normal.
7. **The second amino acid does not appear in the model at all.** Equations 1-3 contain only R, AA
   and D. Every competitor effect is absorbed into refitted k_F and k_E. This is exactly the
   "competition as a multiplier" pattern the registry's policy 6 refuses; the repository is right
   to fit named mass-action channels against these numbers rather than transcribe them.
8. **Amino-acid stock impurities are raised by the authors and not measured here.** They cite
   Becalski et al. for 0.3 % asparagine in an aspartic acid standard and 0.8 % cysteine in a
   glutamine standard, and argue their own glutamine (>=99.5 %) is too pure for contamination to
   explain a 3.2x effect. The argument is sound but the impurity of their own stocks was not
   assayed.
9. **What this paper does NOT contain**: any pH series; any water-activity or moisture series; any
   second laboratory; any absolute level; any sampling-time list; any identification of the
   elimination products (D is written as a three-way "or"); any measurement of asparagine, glucose
   or the competitor amino acid (only acrylamide is quantified); any supplementary material.
10. **What to request from the authors**: (i) the acrylamide concentration-time data behind Figures
    2a-e (five systems x four temperatures x nine times); (ii) the sampling times; (iii) the full
    covariance matrix, given the 0.94-0.98 k_F/k_E correlation; (iv) a copy of Claeys 2005b
    (JAFC 53:9999-10005), which is this scheme's source paper and is **not on disk** — it is the
    obvious next acquisition for this lane.
11. **Registry gaps against `data/keys/compounds.yml`**: `acrylamide` present. **Asparagine,
    glucose, glutamine, cysteine, lysine and alanine absent.** The lane's `Gln`, `Lys`, `Ala`,
    `Cys` species are network-local names, not registry ids; a benchmark or a competitor-panel row
    built from this paper would need all six keyed.
