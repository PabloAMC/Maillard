# Knol 2005 — EXTRACTION (equimolar glucose + asparagine 0.2 mol/L in 0.1 M phosphate pH 6.8, aqueous, 120/140/160/180/200 C, 0-45 min; six-step multiresponse network with an explicit acrylamide DEGRADATION step; 30 rate constants + 6 activation energies with 95 % HPD)
### A SECOND LABORATORY on the acrylamide trunk: Wageningen's aqueous glucose-asparagine network prints the condensation, the acrylamide-forming step AND the first-order acrylamide elimination at five temperatures with intervals, all at the same T_ref = 160 C the shipped lane uses — and the elimination constant, 88.1e-3 min^-1, is the one number in the corpus that lets Leuven's elimination be checked by somebody else.

**Source on disk:** `data/articles/knol2005.pdf` (7 pp., J. Agric. Food Chem. 53:6133-6139).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/knol2005.txt`); Tables 1 and 2 came through clean and are re-typed in full
below. Figures 1A/1B (the reaction networks), 2A-E (experimental concentration-time data at five
temperatures) and 3A-E (model fit over the same data) are images: **every concentration-time datum
in this paper is figure-only.** There is no supplementary material. Repo status before this
dossier: Knol 2005 is cited by `src/kinetic_core/parameters_acrylamide.py` (module docstring,
`REFUSED_PARAMETERS`, and the fitted-step conditions string) and by
`k3_final_parameter_inventory.md` §A.2, but has **no extraction dossier** — its numbers reach the
repository only through that inventory's summary rows.

## 0. Identity

| field | value |
|---|---|
| Title | "Toward a Kinetic Model for Acrylamide Formation in a Glucose-Asparagine Reaction System" |
| Authors | Jeroen J. Knol, Wil A. M. van Loon (these two contributed equally), Jozef P. H. Linssen (corresponding), Anne-Laure Ruck, Martinus A. J. S. van Boekel, Alphons G. J. Voragen — Product Design and Quality Management Group **and** Laboratory of Food Chemistry, Wageningen University, The Netherlands |
| Venue | J. Agric. Food Chem. 2005, 53 (15), 6133-6139. Received 7 March 2005, revised 20 May 2005, accepted 21 May 2005, web 28 June 2005 |
| DOI / article ID | 10.1021/jf050504m (printed as `JF050504M`) |
| Naming | "Decarb. Amadori comp." = decarboxylated Amadori compound; "Product X" = "assumed product(s) formed from acrylamide" (never identified); melanoidins reported as **moles of sugar incorporated** via Lambert-Beer at 470 nm with eps = 282 L mol^-1 cm^-1 |
| Lineage | network taken from Stadler et al. 2002 / 2004 (the N-glycoside route, NOT the Mottram-Wedzicha Strecker route); multiresponse machinery from Martins & van Boekel 2005 and van Boekel 1996; eps from Leong 1999 (Leeds PhD) |
| Companions on disk | `knol2009_extraction.md` (potato crisps, empirical), `knol2010_extraction.md` (fructose, pH 5.5), `claeys2005_extraction.md` and `devleeschouwer2009_extraction.md` (the Leuven lane the repository is fitted on), `martins2005_extraction.md` (the glucose-glycine parent network) |

## 1. Why it matters

The shipped acrylamide lane is wave B3. Its constants live in
`src/kinetic_core/parameters_acrylamide.py` and its network in `src/kinetic_core/acrylamide.py`.
Every populated `MEASURED_ACRYLAMIDE` row in that registry comes from **one laboratory** —
Leuven (Hendrickx): Claeys 2005 for the dilute-aqueous lumps, De Vleeschouwer 2009 Part I
(`_DV1_SOURCE`) for the glucose trunk, De Vleeschouwer 2009 Part II (`_DV2_SOURCE`) for the
cysteine channels. The scorecard's stated gap is "a second laboratory's constants; real-food
matrices" (`results/validation/data_wishlist.md` and the path-by-path scorecard figure).

**This paper is that second laboratory, on the same three steps.** It fits, by multiresponse
regression on five temperatures simultaneously, exactly the slots the repository's network carries:

| repository step (`acrylamide.py`) | registry key | this paper's step |
|---|---|---|
| `a_asn_glc_sb`: Asn + Glc -> SBA | `k_asn_glc` (order 2) | **k1**, glucose + asparagine -> Schiff base |
| (the same slot, fructose) | — | **k3**, fructose + asparagine -> Schiff base |
| `a_int1_acr`: INT1 -> ACR | `k_int1_acr` (order 1) | **k4**, Schiff base / decarboxylated Amadori -> acrylamide |
| `a_int1_mel`: INT1 -> melanoidin | `k_int1_mel` (**fitted**, no literature value accepted) | **k5**, Schiff base -> melanoidins |
| `a_acr_dp`: ACR -> degradation products | `k_acr_dp` (**fitted** across three labs) | **k6**, acrylamide -> Product X, first order |

The ELIMINATION step is the one the repository most needs an outside opinion on. The lane's own
note (policy 5 in `parameters_acrylamide.py`) records that the retired FAST lane had **no**
elimination at all and was ~40x under-responsive to time, and the panel currently has acrylamide's
time shape inverted against a measurement (a 0.5 % glucose + 0.5 % asparagine solution measured at
28, 912 and 1459 ppb after 10, 20 and 30 min, where the core predicts a fall;
`results/legacy_lane/validation/maillard_path_holdout_frozen_predictions.json`). **k6 is a
first-order acrylamide elimination constant at five temperatures, each with a 95 % HPD interval, in
an aqueous pot at pH 6.8.** It is the only such series in the corpus from outside Leuven.

Two structural findings also bear on the network as built. (a) The authors' model discrimination
**rejects** the reverse isomerisation fructose -> glucose (Table 2), the same simplification De
Vleeschouwer 2009 Part I later adopts. (b) There is no Amadori/Schiff resolution here either: the
"Schiff base" is a single unmeasured lump feeding both acrylamide and melanoidins, which is exactly
the topology `ASN_SCHIFF_AMADORI_SPLIT` declares as pinned-and-inert.

What this paper does NOT give the repository: any concentration-time number (all figure-only), any
pH variation (one pH), any water-activity variation (aqueous only), any identification of the
acrylamide degradation products, and — printed on the table header — the correct unit for its two
second-order constants (section 3, Flags 1).

## 2. Methods as they matter to a model

- **Pot.** "Equimolar solutions of glucose and asparagine (0.2 M) were prepared in phosphate buffer
  (0.1 M, pH 6.8)." So **D-glucose 200 mmol/L and L-asparagine 200 mmol/L**, phosphate 100 mmol/L,
  initial pH 6.8. Chemicals: D-glucose, D-fructose, Na2HPO4, KH2PO4 (Merck), L-asparagine (Fluka),
  all analytical grade. **Water activity ~1.0** (dilute aqueous) — this is the axis on which the
  paper differs most from the shipped constants, which are a freeze-dried powder at a_w 0.92.
- **Vessel and heating.** "Samples (10 mL) were heated in hermetically closed screw-capped glass
  tubes (Schott, 16 x 160 mm) at 120, 140, 160, 180, and 200 C in an oil bath. The tubes were
  immersed in the oil up to the cap." Sampling at **0, 1, 2, 4, 9, 15, 30 and 45 min**; immediate
  ice cooling; stored at -20 C. **Triplicate.** The headspace (10 mL sample in a ~32 mL tube) is
  air and is not described; no thermocouple trace is reported and **no heat-up correction is
  applied** — the fit treats every run as isothermal from t = 0 (contrast Claeys 2005 and
  De Vleeschouwer 2009, both of which integrate a logged temperature-time profile).
- **pH.** Set to 6.8 by the buffer at t = 0 and **never re-measured**. The companion paper Knol 2010
  measures the drift in the same kind of pot and finds it large (Fig. 2G there). Treat 6.8 as an
  initial pH.
- **Acrylamide.** HPLC-UV after Barber et al. 2001: sample diluted 1:10 in water, centrifuged 5 min
  at 16000g; Synergi 4 um Hydro-RP C18 (80 A, 250 x 2.00 mm) + AJO-4286 guard; 20 uL injected;
  isocratic 1 % methanol / 99 % 5 mM heptanesulfonic acid, 0.2 mL/min, 20 C; **detection by
  absorbance at 200 nm**, t_r = 6.2 min; **external standard**, calibration curve; **LOD
  10 ug/kg**. No internal standard, no mass spectrometry, no isotope dilution (contrast Claeys 2005
  and De Vleeschouwer, both GC-MS with methacrylamide + butyramide internal standards).
- **Sugars.** Same HPLC, ION-300 ion-exchange column at 85 C, 2.5 mM H2SO4 at 0.4 mL/min,
  refractive index, external standard. Mannose was looked for and **not found**; fructose was.
- **Asparagine.** Diluted 1:1000; EZ:faast kit (Phenomenex); norvaline internal standard (20 nmol);
  Carlo Erba GC5300, Zebron amino acid column 10 m x 0.25 mm, injection 250 C split 1:15, oven
  110 -> 250 C at 20 C/min + 1 min hold; external standard.
- **Melanoidins.** A470 on a UV-1601, diluted as needed, Lambert-Beer with **eps = 282 L mol^-1
  cm^-1** (Leong 1999, for **glucose/asparagine** melanoidins), so the quantity is *moles of sugar
  incorporated into brown polymer*. The authors record that at 160, 180 and 200 C insoluble
  particles formed, scattered light, and **caused the absorbance to fall and the melanoidin
  concentration to be under-estimated**; at 45 min / 160-200 C, 30 min / 180-200 C and 15 min /
  200 C no reliable values were obtained.
- **Fitting.** Athena Visual Studio v10.0; differential equations by mass action, numerical
  integration; **non-linear regression with the determinant criterion** on the individual triplicate
  concentrations (not the means); Arrhenius reparameterised as k = X exp(-Y Ea) with
  X = k0 exp(-Ea/(R T_av)) and Y = (1/R)(1/T - 1/T_av), T_av = mean of the five temperatures.
  **T_av = 160 C exactly**, so the printed 160 C column IS the estimated X. All five temperatures
  fitted simultaneously. 95 % highest-posterior-density intervals reported. Model discrimination by
  corrected Akaike criterion and posterior probability (Table 2).
- **Reference temperature.** T_av = 160 C = **433.15 K**, the same `T_REF_A_K` the registry uses.
  No conversion is needed to put this paper's constants next to the shipped ones.

## 3. Tables re-typed

### Table 1. "Estimates of Rate Constants and Activation Energies as Found by Kinetic Modeling for the Proposed Kinetic Model"

Header exactly as printed: **`k (x10^-3 min^-1)`** for all six rows; footnote `a`: "±95 % highest
posterior density (HPD) interval." (On the unit of k1 and k3 see Flags 1 — the header is
incomplete, and the same laboratory's Knol 2010 prints the corresponding constant as
`10^-3 L mmol^-1 min^-1`.)

| k | 120 C | 140 C | 160 C | 180 C | 200 C | Ea (kJ/mol) |
|---|---|---|---|---|---|---|
| k1 | 0.131 ± 0.025 | 0.308 ± 0.049 | 0.668 ± 0.13 | 1.35 ± 0.36 | 2.58 ± 0.89 | 57.6 ± 8.0 |
| k2 | 4.98 ± 1.2 | 16.7 ± 3.0 | 50.1 ± 9.6 | 136 ± 35 | 341 ± 114 | 81.7 ± 8.6 |
| k3 | 0.0819 ± 0.045 | 0.369 ± 0.14 | 1.45 ± 0.44 | 5.04 ± 1.6 | 15.8 ± 6.4 | 102 ± 14 |
| k4 | 0.176 ± 0.081 | 0.712 ± 0.23 | 2.53 ± 0.51 | 8.05 ± 1.2 | 23.2 ± 4.3 | 94.4 ± 11 |
| k5 | 15.7 ± 3.5 | 28.4 ± 4.6 | 48.7 ± 5.9 | 79.6 ± 8.9 | 125 ± 16 | 40.1 ± 5.0 |
| k6 | 7.96 ± 5.1 | 28.1 ± 13 | 88.1 ± 25 | 250 ± 45 | 650 ± 136 | 85.1 ± 14 |

**Which step is which.** Figure 1B, which carries the numbering, is an image. The running text
names four of the six directly and unambiguously:

- **k1** — "the rate constant of reaction 1 (k1)" is compared throughout against k3 as "the same
  formation reaction with glucose and asparagine as reactants": **Glc + Asn -> Schiff base**.
- **k3** — "the rate of the Schiff base formation reaction from fructose and asparagine (k3)":
  **Fru + Asn -> Schiff base**.
- **k6** — "the model was able to fit the loss of acrylamide (k6) to the experimental
  observations": **acrylamide -> Product X**, first order.
- The network is described in words: "the reaction between asparagine and fructose was added ...
  The formation of melanoidins was suggested to result from further reaction of Schiff base, and
  the loss of acrylamide was accounted for by the putative formation of hitherto unknown
  product(s)", plus "the reversible isomerization reaction between glucose and fructose".

That leaves **k2**, **k4** and **k5** to the three remaining steps of Figure 1B: glucose ->
fructose isomerisation, Schiff base -> acrylamide, and Schiff base -> melanoidins. The assignment
below is **mine, by elimination and magnitude**, and is not printed:

- **k2 = glucose -> fructose (isomerisation)** — the only step of the three that must be fast and
  strongly temperature-dependent enough to explain the measured fructose maxima at 9 / 4 / 2 min
  (160 / 180 / 200 C); it is the second largest constant in the table.
- **k4 = Schiff base -> acrylamide** — its 160 C value 2.53e-3 min^-1 sits within a factor 1.4 of
  De Vleeschouwer 2009 Part I's k_Fref (3.57e-3 min^-1, same slot, same T_ref) and within 2 of
  Knol 2010's k3 (5.0e-3 min^-1, the step the 2010 text *names* as "the formation of acrylamide
  from the Schiff base"). Knol 2010's own numbering runs k1 isomerisation, k2 condensation, k3
  Schiff -> acrylamide, k4 acrylamide degradation, k5 Schiff -> melanoidins; if 2005's numbering is
  the same family, k4 and k5 here are the two Schiff-base fates, which is what the magnitudes say.
- **k5 = Schiff base -> melanoidins** — the larger of the two Schiff-base fates (48.7e-3 vs
  2.53e-3 min^-1 at 160 C), which is required: the acrylamide yield in this pot is ~1.5 % of the
  initial asparagine, so the browning fate must take the great majority of whatever condenses.

**k5/k4 = 19.2 at 160 C (mine)**, i.e. the intermediate goes to browning 19 times more often than
to acrylamide — a 4.9 % acrylamide share of the intermediate's flux, against a 1.5 % share of the
initial asparagine, the difference being the asparagine that never condenses. That ratio is what
the repository's `k_int1_mel` is fitted to represent and cannot
take from De Vleeschouwer (whose k_M is marked "NO PHYSICAL MEANING" and is refused). If the
k4/k5 assignment above is right, **this paper prints the partition the repository currently
fits.** It should be treated as a candidate, not a fact, until Figure 1B is read from the page
image.

### Table 2. "Model Discrimination Results for Our Proposed Model with (A) or without (B) the Reversible Isomerization Reaction from Fructose to Glucose"

| model | p | SS | n | AICc | ΔAICc | PPB | PPS |
|---|---:|---|---:|---:|---:|---:|---:|
| A | 16 | 1.32 x 10^5 | 600 | 6505.2 | 27.3 | -89.83 | 0.234 |
| B | 14 | 1.29 x 10^5 | 600 | 6477.9 | 0.0 | -89.33 | 0.757 |

Footnote: "p, number of parameters; SS, residual sum of squares; n, number of data points including
the replicates; AICc, Akaike criterion; ΔAICc, AICc difference with the smallest value taken as
reference; PPB, log10 of posterior probability; PPS, normalized posterior probability share."
Verdict in the text: "The results from our test ... support the model without the reversible
isomerization reaction." **Table 1 therefore reports model B.**

### Numbers printed in the running text (everything else in this paper is figure-only)

| quantity | value | where |
|---|---|---|
| complete loss of glucose | at 9 min (180 C), 5 min (200 C) | Results, Fig. 2A |
| complete loss of asparagine | at 30 min (180 C), 9 min (200 C) | Results, Fig. 2B |
| fructose maximum | at 9 / 4 / 2 min (160 / 180 / 200 C); complete loss after 30 / 9 min at 180 / 200 C | Results, Fig. 2C |
| mannose | not formed (looked for, absent) | Results |
| acrylamide time shape | rises then steady state at 140 and 160 C; rises then "fast decrease" at 180 and 200 C; the maximum coincides with sugars reaching zero | Results, Fig. 2D |
| k3 vs k1 | k3 is 35 % lower than k1 at 120 C; 2-3x higher than k1 at 180 C | Results |
| melanoidin absorbance invalid | 45 min at 160/180/200 C; 30 min at 180/200 C; 15 min at 200 C | Results |
| average precision | "about ±27 % for the rate constants and ±11 % for the activation energies" | Results |
| eps (melanoidin, Glc/Asn) | 282 L mol^-1 cm^-1 | Methods (from Leong 1999) |
| acrylamide LOD | 10 ug/kg | Methods |
| acrylamide yield in this system | **not printed here**; Knol 2010's conclusion quotes it as **1.5 % of the initial asparagine** at pH 6.8 with glucose | Knol 2010 conclusion |

**Concentration-time data: FIGURE-ONLY.** Figures 2A-E (glucose, asparagine, fructose, acrylamide,
melanoidins; five temperatures; means ± SD of triplicates) and 3A-E (the same with fitted lines) are
the only place any concentration appears. Per house rule they are not typed as numbers.

### Arithmetic on the printed constants (all mine)

**1. Ratios inside the table at 160 C (unit-independent within a column of the same order).**
k6/k4 = 88.1/2.53 = **34.8** — the acrylamide sink is ~35x its own source constant, so acrylamide
in this pot is a fast-turning-over intermediate whose level is set by (k4 x [Schiff]) / k6, not by
k4 alone. k5/k4 = **19.2** (the browning-vs-acrylamide partition, subject to the k4/k5 assignment
above). k3/k1 at 160 C = 1.45/0.668 = **2.17**, i.e. fructose is the more reactive partner at
T_av, rising to 3.7 at 200 C and falling to 0.63 at 120 C (the paper's own "35 % lower at 120 C").

**2. Half-life of acrylamide from k6 (mine, first order).** t_1/2 = ln2/k6: **87 min at 120 C,
24.7 min at 140 C, 7.9 min at 160 C, 2.8 min at 180 C, 1.07 min at 200 C.** This is the single most
useful derived number in the paper for the panel's inverted-time-shape problem: at 120 C an
acrylamide pool loses only ~21 % in 30 min, so a *measured* rise over 10-30 min is fully
compatible with this laboratory's elimination constant, whereas at 160 C the same pool would be
past its maximum well before 30 min.

**3. Cross-check of the Arrhenius fit (mine).** From k6 at 120 and 200 C:
Ea = R ln(650/7.96) / (1/393.15 - 1/473.15) = 8.314e-3 x 4.400 / 4.302e-4 = **85.0 kJ/mol**,
against the printed 85.1 ± 14. The same two-point check reproduces k1 (57.4 vs 57.6), k2 (81.7 vs
81.7), k3 (101.9 vs 102), k4 (94.6 vs 94.4) and k5 (40.1 vs 40.1). **The table is internally
consistent to better than 0.5 %**, which is a strong indication that the five temperature columns
are derived from X and Ea rather than fitted independently — and therefore that the five columns
carry only two degrees of freedom each, not five.

**4. The 129 kJ/mol barrier is not here.** The largest activation energy in this paper is
**102 ± 14 kJ/mol** (k3), whose upper 95 % bound is 116. The acrylamide-forming k4 is
**94.4 ± 11** and the elimination k6 is **85.1 ± 14**. `REFUSED_PARAMETERS` in
`parameters_acrylamide.py` states this and it is confirmed here from the primary table: **no
number in Knol 2005 is 129, and none can round to it.** The registry's paired claim that
`safety_reference_payloads.json` entries[27] mis-states the pair as 52.1 / 72.9 where the true
values are 94.4 ± 11 and 85.1 ± 14 is also confirmed against the printed table.

**5. Parameter count does not close (mine).** Table 1 prints 6 rate constants and 6 activation
energies = 12 parameters, but Table 2 counts model B at p = 14 and model A at p = 16. The
difference A - B = 2 is exactly one (k, Ea) pair, consistent with removing one reversible step. The
residual 2 parameters in model B are **not identified anywhere in the text or the tables**. A
plausible reading is two fitted initial concentrations, but that is a guess and is recorded as
Flag 4, not as a number.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** `acrylamide` is keyed (id `acrylamide`, CAS
79-06-1). **Asparagine, glucose, fructose, aspartic acid, melanoidins and "Product X" are NOT in
the registry** — the registry is a product/marker list and carries no Maillard reactants (checked:
75 ids, none of them a sugar or an amino acid). Every row below shares: 200 mmol/L glucose +
200 mmol/L asparagine, 0.1 M phosphate, **initial** pH 6.8, dilute aqueous (a_w ~ 1.0), sealed
glass tube with air headspace, unstirred, 120-200 C, 0-45 min, triplicate, fitted with no heat-up
correction, T_ref = 160 C = 433.15 K.

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| Glc + Asn -> Schiff base | k1 at 120/140/160/180/200 C | 0.131 ± 0.025 / 0.308 ± 0.049 / **0.668 ± 0.13** / 1.35 ± 0.36 / 2.58 ± 0.89 | printed `x10^-3 min^-1`; **read as 10^-3 L mmol^-1 min^-1 = M^-1 min^-1** (Flags 1) | aqueous pH 6.8, a_w ~1.0 | **second order**, mass action in [Asn][Glc] | Table 1 p. 6137 | measured_rate (unit corrected by inference) |
| " | Ea(k1) | 57.6 ± 8.0 | kJ/mol | 120-200 C | — | Table 1 | measured_barrier |
| Glc -> Fru (isomerisation) | k2 at the five T | 4.98 ± 1.2 / 16.7 ± 3.0 / **50.1 ± 9.6** / 136 ± 35 / 341 ± 114 | 10^-3 min^-1 | as above | first order in glucose | Table 1 | measured_rate (step identity mine, Flags 3) |
| " | Ea(k2) | 81.7 ± 8.6 | kJ/mol | " | — | Table 1 | measured_barrier |
| Fru + Asn -> Schiff base | k3 at the five T | 0.0819 ± 0.045 / 0.369 ± 0.14 / **1.45 ± 0.44** / 5.04 ± 1.6 / 15.8 ± 6.4 | printed `x10^-3 min^-1`; **read as 10^-3 L mmol^-1 min^-1 = M^-1 min^-1** | as above | **second order** | Table 1 | measured_rate (unit corrected by inference) |
| " | Ea(k3) | 102 ± 14 | kJ/mol | " | — | Table 1 | measured_barrier (**the largest in the paper**) |
| Schiff base -> acrylamide | k4 at the five T | 0.176 ± 0.081 / 0.712 ± 0.23 / **2.53 ± 0.51** / 8.05 ± 1.2 / 23.2 ± 4.3 | 10^-3 min^-1 | as above | first order in the unmeasured intermediate | Table 1 | measured_rate (step identity mine, Flags 3) |
| " | Ea(k4) | 94.4 ± 11 | kJ/mol | " | — | Table 1 | measured_barrier |
| Schiff base -> melanoidins | k5 at the five T | 15.7 ± 3.5 / 28.4 ± 4.6 / **48.7 ± 5.9** / 79.6 ± 8.9 / 125 ± 16 | 10^-3 min^-1 | as above | first order | Table 1 | measured_rate (step identity mine, Flags 3) |
| " | Ea(k5) | 40.1 ± 5.0 | kJ/mol | " | — | Table 1 | measured_barrier |
| **acrylamide -> Product X (ELIMINATION)** | **k6 at the five T** | 7.96 ± 5.1 / 28.1 ± 13 / **88.1 ± 25** / 250 ± 45 / 650 ± 136 | **10^-3 min^-1** | as above | **first order in acrylamide** | Table 1 | measured_rate — **carry the authors' caveat verbatim: "the model was not restrained by experimental data for the products formed in the degradation reaction"** |
| " | Ea(k6) | **85.1 ± 14** | kJ/mol | 120-200 C | — | Table 1 | measured_barrier |
| acrylamide half-life from k6 | t_1/2 | 87 / 24.7 / 7.9 / 2.8 / 1.07 | min at 120/140/160/180/200 C | as above | — | derived from Table 1 (mine) | derived_assumption (arithmetic only; assumes first order and no re-formation) |
| browning-vs-acrylamide partition | k5/k4 | 89.2 / 39.9 / **19.2** / 9.9 / 5.4 at 120/140/160/180/200 C | — | as above | — | derived from Table 1 (mine) | within_study_ratio (rests on the k4/k5 assignment, Flags 3) |
| sink-vs-source on acrylamide | k6/k4 | 45 / 39 / **34.8** / 31 / 28 at 120-200 C | — | as above | — | derived (mine) | within_study_ratio |
| fructose-vs-glucose reactivity | k3/k1 | 0.63 / 1.20 / **2.17** / 3.73 / 6.12 at 120-200 C | — | as above | — | derived (mine) | within_study_ratio |
| eps, melanoidins from Glc/Asn | 282 | L mol^-1 cm^-1 | A470, aqueous | — | Methods (from Leong 1999) | level_only (**not measured here**; borrowed) |
| model choice | reverse isomerisation Fru -> Glc | rejected (ΔAICc 27.3; PPS 0.757 vs 0.234) | — | — | — | Table 2 | within_study_ratio (structural) |
| acrylamide, glucose, fructose, asparagine, melanoidin concentration-time courses at five temperatures | — | — | mmol/L vs min | — | — | Figs. 2A-E, 3A-E | **figure_only** |

### Can these be put on the same basis as the trunk's constants? Step by step.

The registry's own T_ref is 160 C and so is this paper's T_av, so **no temperature transport is
needed** — the 160 C column is directly comparable to `k_ref`. What differs is the matrix
(dilute aqueous 0.2 mol/L, a_w ~1.0 here vs freeze-dried powder at a_w 0.92 and ~2.7 mol/kg there)
and, for k1/k3, the printed unit.

**(a) The condensation, second order — comparable, and it agrees.** Reading k1 as
0.668e-3 L mmol^-1 min^-1 = **0.668 M^-1 min^-1** (Flags 1), the shipped `k_asn_glc` is
De Vleeschouwer 2009 Part I's k_INTg = **1.70 ± 1.05 M^-1 min^-1** at the same 160 C. Knol's value
is **2.5x lower** and sits **inside** De Vleeschouwer's 95 % HPD (lower bound 0.65). A second-order
constant is the one kind that transports across a concentration change, so this is a genuine
cross-laboratory, cross-matrix agreement on the initiation step — the single most valuable line in
this dossier for the scorecard's stated gap. **The barriers do not agree**: 57.6 ± 8.0 (here)
against 117.5 ± 25.2 (Part I), intervals separated by 27 kJ/mol.

**(b) The acrylamide-forming step, first order — comparable in unit, but the lump differs.**
k4 = 2.53e-3 min^-1 at 160 C against the shipped `k_int1_acr` = 3.57e-3 min^-1 (Part I k_Fref,
whose printed unit `10^-3 mm^-1` the inventory already identifies as a typo for min^-1). **A
factor of 1.4.** Both are first order in an unmeasured intermediate, but "Schiff base" here and
"Int1" there are different lumps in different networks, so the agreement is suggestive rather than
decisive. Barriers: 94.4 ± 11 here against 159.2 ± 29.5 there — **not overlapping** (upper 105.4 vs
lower 129.7).

**(c) The elimination, first order — fully comparable, and this is the high-value row.**
k6 = **88.1e-3 min^-1 at 160 C, Ea 85.1 ± 14 kJ/mol**. The two Leuven values at the same T_ref are
Claeys 2005 control k_Eref = **111.1e-3 min^-1, Ea 167.21 ± 4.30** and De Vleeschouwer 2009 Part I
glucose k_Eref = **0.10 min^-1, Ea 113.2 ± 32.3**. The three rates span **1.26x** (0.0881 to
0.1111); the three barriers span **2.0x** and Knol's interval does not overlap Claeys' (99.1 vs
162.9) though it does overlap De Vleeschouwer's (99.1 vs 80.9). This reproduces exactly the picture
`parameters_acrylamide.py` records in its GROUP 2 comment — and it now rests on the primary
tables rather than on the inventory's summary. **A re-calibration wave can fit one (k_ref, Ea) pair
to all three and let the residual show the barrier conflict; it cannot pick a barrier from this
paper alone.**

**(d) What CANNOT be transported.** Nothing here measures water activity, pH, a competing amino
acid, or a real food. The a_w gap between this pot (~1.0) and the registry's declared
`aw_of_measurement` (0.92) is 0.08 — inside the 0.1 warning band the registry emits — but that is
an accident of the numbers, not evidence: the a_w windows in `acrylamide_conditions.py` are
De Vleeschouwer's, and this paper adds no point to them. Every concentration is figure-only, so
this paper supplies **no benchmark row**, only constants.

## 5. Flags

1. **The unit on Table 1's header is wrong for k1 and k3, and the correction is licensed by the
   same laboratory.** The header reads `k (x10^-3 min^-1)` for all six rows, but k1 and k3 are
   bimolecular steps (sugar + asparagine) fitted by mass action, which cannot carry min^-1. Knol
   2010, same group, same first author, same network family, prints its condensation constant with
   an explicit second footnote: `** (10^-3 l mmol^-1 min^-1)` against `* (10^-3 min^-1)` for
   everything else. The 2010 text then compares the two directly — "The rate constant for the
   formation of the Schiff base from fructose and asparagine (k2) at pH 5.5 is lower than that at
   pH 6.8 (Knol et al., 2005)" — which only parses if the two are in the same unit (2010 k2 at
   160 C = 0.27e-3 vs 2005 k3 = 1.45e-3, a factor 5.4 lower, as claimed). **Working reading:
   k1 and k3 are in 10^-3 L mmol^-1 min^-1, i.e. numerically equal to M^-1 min^-1.** This is an
   inference, marked as such in section 4, and it is the difference between agreeing with the
   shipped `k_asn_glc` within its interval and being 2500x below it.
2. **The acrylamide degradation is unconstrained by data and the authors say so.** "For the
   estimation of the acrylamide concentration, the model was not restrained by experimental data
   for the products formed in the degradation reaction, and therefore, the model was able to fit
   the loss of acrylamide (k6) to the experimental observations. More knowledge about the reaction
   products from acrylamide would validate the estimated rate constants for the formation and loss
   of acrylamide." k6's own HPD at 120 C is ±64 % (7.96 ± 5.1). The constant is a *measurement of
   how fast acrylamide disappears*, not of a named reaction; carry it as `k_acr_dp`'s evidence,
   never as a mechanism.
3. **k2, k4 and k5 are assigned to steps by my inference, not by this paper — but Knol 2010
   corroborates two of the three.** Figure 1B, the only place the numbering is drawn, is an image;
   the text names only k1, k3 and k6 and describes the network in prose. The assignment in
   section 3 (k2 isomerisation, k4 Schiff -> acrylamide, k5 Schiff -> melanoidins) rests on
   magnitude and on the prose. Two independent checks from Knol 2010 (same group, same first
   author) support it: (i) Knol 2010 writes that "the activation energy for the formation of
   acrylamide from the Schiff base was almost a factor of 3 lower" than at pH 6.8, and its own
   Ea(k3) = 35 against **this paper's Ea(k4) = 94.4** gives 2.7 — **so k4 here is the
   acrylamide-forming step**; (ii) Knol 2010's Schiff -> melanoidin step carries Ea = 40 ± 26
   against **this paper's Ea(k5) = 40.1 ± 5.0**, the same value to one decimal. The remaining
   assignment, k2 = isomerisation, is by elimination only. **Read Figure 1B from the page image
   before any of k2, k4 or k5 enters a registry.** k6 and the two condensation constants do not
   depend on this. Note that Knol 2010's *rate* comparison of the same step ("about a factor of 2
   higher ... at 120 C") does **not** reproduce against the printed values: 1.9e-3 there against
   0.176e-3 here is 10.8x (mine).
4. **The parameter count does not close.** Table 1 shows 12 estimated parameters; Table 2 says
   model B has 14 and model A 16. Two parameters of model B are unaccounted for in the paper.
   Ask the authors, or treat the reported HPD intervals as conditional on two unnamed nuisance
   parameters.
5. **No heat-up correction.** The tubes were immersed in an oil bath and the fit is isothermal from
   t = 0. Claeys 2005 and both De Vleeschouwer 2009 papers integrate a 2-4 s thermocouple log
   through the heat-up and cool-down. In a 10 mL glass tube the heat-up to 200 C is minutes, and
   the earliest sample is at 1 min. **The 180 and 200 C columns are the ones this bites**, and it
   biases those k downward, which flattens the fitted Arrhenius slope — a plausible partial
   explanation for why every barrier here is lower than the corresponding Leuven barrier.
6. **The five temperature columns are not five independent measurements.** My two-point Arrhenius
   check reproduces all six printed Ea to better than 0.5 %, so each row carries two free
   parameters (X and Ea) presented as six numbers. Do not treat the 120 C and 200 C columns as
   additional evidence when weighting a fit.
7. **Acrylamide by HPLC-UV at 200 nm with an external standard.** No internal standard, no MS, no
   isotope dilution; a 200 nm detection in a browning aqueous Maillard pot is the least selective
   acrylamide method in this five-paper set (Claeys and De Vleeschouwer both use GC-MS/PCI with
   methacrylamide and butyramide internal standards). Systematic level error is plausible; the
   *shape* the kinetics is fitted to is more robust than the absolute level.
8. **pH is initial only and never re-measured**, in an unbuffered-in-practice pot (0.1 M phosphate
   against 0.2 M reactants). Knol 2010 measures the same drift and finds it goes down by more than
   a unit and then partly back up. `acrylamide_conditions.py` installs pH factors indexed to
   *initial* pH for exactly this reason; this paper's 6.8 must be read the same way.
9. **The melanoidin response is under-estimated at 160-200 C by the authors' own account** (light
   scattering by insoluble particles), and its eps = 282 L mol^-1 cm^-1 is borrowed from a Leeds
   PhD thesis, not measured here. Since k5 is fitted against that response, **k5 inherits the
   under-estimation**, and so does the k5/k4 partition derived from it.
10. **What this paper does not contain**: any pH series; any water-activity or moisture series; any
    competing amino acid; any real-food matrix; any identification of Product X; any tabulated
    concentration; any acrylamide yield figure (the 1.5 % figure comes from Knol 2010's
    conclusion); any supplementary material.
11. **What to request from the authors**: (i) the numeric concentration-time data behind Figures
    2A-E (five temperatures x eight times x five responses x three replicates); (ii) Figure 1B at
    resolution, or simply the k-to-step mapping; (iii) confirmation of the unit on k1 and k3;
    (iv) the identity of the two extra parameters in Table 2's p = 14.
12. **Registry gaps against `data/keys/compounds.yml`**: `acrylamide` is present. **Asparagine,
    glucose, fructose, the Schiff base / decarboxylated Amadori lump, melanoidins and Product X
    are all absent.** The lane's internal species names (`Asn`, `Glc`, `SBA`, `INT1`, `ACR`,
    `MEL_C`, `MEL_N`) are network-local and are not registry ids; a benchmark built from this paper
    would need at least asparagine and glucose keyed.
