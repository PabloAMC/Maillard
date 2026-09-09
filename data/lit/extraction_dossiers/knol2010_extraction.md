# Knol 2010 — EXTRACTION (equimolar fructose + asparagine 0.1 mol/L in 0.1 M phosphate initial pH 5.5, aqueous, 120/140/160/180/200 C, 0-64 min; seven-step multiresponse network with acetic and formic acid; 35 rate constants + 7 activation energies with 95 % HPD — and NO acrylamide elimination step, because its barrier came out negative and the step was deleted)
### The second laboratory's condensation constant printed in a unit the repository can use without inference — 0.27 x 10^-3 L mmol^-1 min^-1 at 160 C for fructose + asparagine, against De Vleeschouwer's 0.22 M^-1 min^-1 for the same step, a 1.2x agreement across two laboratories and two matrices; and the paper that proves the acrylamide elimination is unidentifiable in this system by deleting it.

**Source on disk:** `data/articles/knol2010.pdf` (11 pp., Food Chemistry 120 (2010) 1047-1057).
Read from the `pdftotext -layout` text layer (`scratchpad/articles/knol2010.txt`), with page 9
re-extracted in `-raw` mode to confirm Table 2's footnoted units. Tables 1 and 2 came through clean
and are re-typed in full below, **including the two unit footnotes, which are the most valuable
thing in the paper.** Figures 1 (temperature profiles), 2A-H (experimental data: fructose,
asparagine, glucose, acrylamide, formic acid, acetic acid, **pH**, melanoidins at five
temperatures), 3A-E (mass balances), 4 (the full Stadler network), 5A-G (model fit) and Schemes 1-4
(the four successive kinetic models) are images: **every concentration-time datum and the entire pH
trajectory are figure-only, and the k-to-step mapping is drawn only in Schemes 3 and 4.** No
supplementary material. Repo status before this dossier: Knol 2010 is named in
`parameters_acrylamide.py` (`REFUSED_PARAMETERS`, the 129 kJ/mol entry and the degradation entry)
and is declared **HOLD-OUT** in `k3_final_parameter_inventory.md` D.5 — "a *third* lab on the same
trunk; the only genuine cross-lab extrapolation test the acrylamide module can have", with a
recommended split (hold out the acrylamide steps, fit the organic-acid and isomerisation steps).
It has **no extraction dossier**.

## 0. Identity

| field | value |
|---|---|
| Title | "Unravelling the kinetics of the formation of acrylamide in the Maillard reaction of fructose and asparagine by multiresponse modelling" |
| Authors | Jeroen J. Knol, Jozef P. H. Linssen (corresponding), Martinus A. J. S. van Boekel — Product Design and Quality Management Group, Wageningen University, The Netherlands |
| Venue | Food Chemistry 120 (2010) 1047-1057. Received 15 September 2008, revised 17 November 2009, accepted 20 November 2009 |
| DOI | 10.1016/j.foodchem.2009.11.049 |
| Naming | "X1" = the (deleted) acrylamide degradation products; "X2" = unidentified asparagine breakdown products; the **Heyns** product is explicitly bypassed ("we bypass the formation of the Heyns product altogether", following Brands & van Boekel 2002a); melanoidins as moles of sugar incorporated, eps = 282 L mol^-1 cm^-1 (Leong 1999, a **glucose/asparagine** value used on a fructose/asparagine pot) |
| Lineage | the Stadler et al. 2004 network (Fig. 4), simplified through four successive schemes; Athena Visual Studio v11.0; the reparameterised Arrhenius of van Boekel 1996; direct descendant of Knol 2005 |
| Companions on disk | `knol2005_extraction.md` (glucose, pH 6.8 — the paper this one is compared against throughout), `knol2009_extraction.md`, `claeys2005_extraction.md`, `devleeschouwer2009_extraction.md` |

## 1. Why it matters

The shipped acrylamide lane (wave B3; `src/kinetic_core/parameters_acrylamide.py`,
`src/kinetic_core/acrylamide.py`) has every one of its populated measured rows from a single
laboratory, Leuven. The scorecard names the gap as "a second laboratory's constants; real-food
matrices" (`results/validation/data_wishlist.md`). Wageningen is the second laboratory, and **this
paper is the only place in the corpus where the second laboratory prints a bimolecular
Maillard-initiation constant in an unambiguous, footnoted second-order unit.**

That matters because of what the registry's `k_asn_glc` is. `parameters_acrylamide.py` calls it
"THE reason this module exists in mass-action form ... the only genuinely second-order Maillard-
initiation constant anywhere in the corpus, which is what lets the network respond to precursor
CONCENTRATION instead of carrying a fixed yield". It is De Vleeschouwer 2009 Part I's k_INTg =
1.70 M^-1 min^-1 (converted once by `per_molar_to_per_mmol` to 1.70e-3 L/(mmol*min)). Nothing has
ever checked it. **Table 2 of this paper prints k2, the fructose + asparagine condensation, as
0.27 ± 0.06 in units of 10^-3 L mmol^-1 min^-1 at exactly the registry's T_ref of 160 C** — i.e.
0.27 M^-1 min^-1 — against De Vleeschouwer 2009 Part I's fructose column k_INTf = **0.22 ± 0.38
M^-1 min^-1**. Two laboratories, two matrices separated by a factor of ~25 in concentration and a
water activity of 1.0 against 0.92, agreeing on the same step within 25 %. Section 4 does that
arithmetic.

The footnote also settles Knol 2005. That paper heads its whole table `k (x10^-3 min^-1)` although
two of its six steps are bimolecular; this paper prints the corresponding step with a **separate
footnote `** (10^-3 l mmol^-1 min^-1)`** and then compares the two constants directly in the text
("The rate constant for the formation of the Schiff base from fructose and asparagine (k2) at
pH 5.5 is lower than that at pH 6.8 (Knol et al., 2005)"). That comparison only parses if the two
are in the same unit. **This paper is the licence for reading Knol 2005's k1 and k3 as
second-order constants** (see `knol2005_extraction.md`, Flags 1).

On **acrylamide ELIMINATION**, this paper's contribution is a negative result and it is a decisive
one. Scheme 3 contained an acrylamide degradation step; "the activation energy for the breakdown
of acrylamide was negative, which suggests that, with increasing temperature, the degradation of
acrylamide would decrease and this is of course not true. Therefore, the degradation route of
acrylamide was omitted in Scheme 4." The authors' explanation: "The degradation of acrylamide in
the fructose-asparagine reaction system (shown in Fig. 2D) was not as significant as that reported
earlier in a glucose-asparagine system with pH 6.8 (Knol et al., 2005) and apparently the model has
difficulty in estimating the breakdown process." **There is therefore NO elimination constant in
this paper**, and `parameters_acrylamide.py`'s `REFUSED_PARAMETERS` entry ("Knol 2010: the Ea went
negative and was deleted") is confirmed verbatim from the primary text. The lane's decision to FIT
`k_acr_dp` rather than transcribe any single source's value is supported by this paper, not
contradicted by it.

Finally, three constants here belong to the trunk rather than to the acrylamide lane and are the
ones the inventory recommends fitting: **the isomerisation (61 ± 8 kJ/mol), the acetic-acid route
(75 ± 10) and the formic-acid route (84 ± 14)**, all with five-temperature rate series. They are
re-typed below with their step identities argued rather than assumed.

## 2. Methods as they matter to a model

- **Pot.** "Equimolar solutions of fructose and asparagine (0.1 M) were prepared in phosphate buffer
  (0.1 M, pH 5.5)." So **D-fructose 100 mmol/L, L-asparagine 100 mmol/L, phosphate 100 mmol/L,
  initial pH 5.5**. Water activity ~1.0. Chemicals from Merck and Fluka, analytical grade.
  Deliberately chosen: "The pH of the system was adjusted to pH 5.5, which is more in line with the
  pH that is found in potatoes (Burton, 1989)." Half the concentration of Knol 2005, and a
  different sugar and pH — the three differences travel together and cannot be separated.
- **pH is NOT constant and the authors say so plainly.** "the buffering capacity of the 0.1 M
  phosphate buffer was not strong enough; dissolving the reactants lowered the initial pH of the
  buffer already from 5.5 to 5.3." The pH then falls further during heating (maximum decrease after
  32 and 16 min at 180 and 200 C), and at 200 C "the final pH even exceeded the initial pH". The pH
  was measured **after cooling to room temperature**, so "the exact pH at the actual temperatures is
  unknown". "For modelling purposes, the change in pH has not been taken into account, which
  implies that possible pH effects, if any, are hidden in the parameters obtained from the models."
  **Every constant in Table 2 is a pH-5.5-initial, pH-drifting constant.** This is the same
  construction `acrylamide_conditions.py` uses for its declared pH factors (indexed to initial pH
  while the pot drifts), and the same caveat applies.
- **Vessel and heating.** "Samples (10 ml) were heated in hermetically closed screw-capped
  **stainless steel** tubes (Workshop, Wageningen University) at 120, 140, 160, 180 and 200 C in a
  heating block (Liebisch)." Sampling at **0, 1, 2, 4, 8, 16, 32 and 64 min**; ice cooling; -20 C
  storage. **Duplicate at 0, 16, 32 and 64 min; triplicate at 1, 2, 4 and 8 min.** Knol 2005 used
  glass tubes in an oil bath and went to 45 min; this one goes to 64 min in steel in a block.
- **Heat-up: measured but NOT corrected.** Fig. 1 shows the temperature profile of the system
  during heating of the steel tubes at all five set-points. "Although the experimental conditions
  were non-isothermal in the first minutes of the experiments (Fig. 1), the system was considered
  isothermal for modelling purposes ... The complexity of the modelling is reduced in this way and
  the temperature dependence is then incorporated in the parameters of the model." **The heat-up is
  therefore absorbed into the constants**, exactly as in Knol 2005 and unlike Claeys 2005 /
  De Vleeschouwer 2009, both of which integrate the logged profile. Fig. 1 is an image; the heat-up
  time is not printed.
- **Acrylamide.** HPLC-UV, adapted from Knol 2005 / Barber 2001. Dilution 1:5 or 1:10, 0.2 um
  filter; Synergi 4 um Hydro-RP C18 (80 A, 250 x 2.00 mm) + AJO-4286 guard; 20 uL; isocratic
  **6 % acetonitrile / 94 % water**, 0.2 mL/min, 20 C; **absorbance at 210 nm**, t_r = 4.8 min;
  **external standard**. (Knol 2005 used 1 % methanol / heptanesulfonic acid and 200 nm.) No
  internal standard, no MS. No LOD printed.
- **Asparagine.** HPLC-**ELSD** (Alltech 3300; drift tube 60 C, gas 2.0 L/min, gain 1); Prevail C18
  5 um 250 x 4.6 mm; gradient of 5 mM heptafluorobutyric acid + 0.7 % TFA against acetonitrile
  (0 % B to 6 min, 0->15 % B from 6 to 8 min, 15->35 % B from 8 to 25 min), 1 mL/min, 20 C;
  external standard. (Knol 2005 used the EZ:faast GC kit with norvaline; the method changed.)
- **Sugars.** HPLC-RI (Martins 2003a), IOA-1000 organic-acids column 300 x 7.8 mm at 85 C, 2.5 mM
  H2SO4 at 0.4 mL/min, external standard. **Mannose was looked for and not found; glucose was
  formed** from fructose.
- **Organic acids.** HPLC-DAD, Prevail organic acids 5 um 250 x 4.6 mm at 20 C, 25 mM KH2PO4
  acidified to pH 2.5 with phosphoric acid, 1.0 mL/min; UV spectrum 180-380 nm for the profile and
  **210 nm for quantification**, external standard.
- **Melanoidins.** A470 on a Cary 50-Bio, Lambert-Beer with **eps = 282 L mol^-1 cm^-1**, the
  Leong 1999 value for **glucose/asparagine** melanoidins applied here to a fructose/asparagine pot;
  the authors note this "was not explored at the high temperatures used in this study" and that
  insoluble high-molecular-weight melanoidins scatter and cause under-estimation at 180 C beyond
  16 min and 200 C beyond 8 min.
- **Fitting.** Athena Visual Studio v11.0; mass-action differential equations, numerical
  integration; non-linear regression by the **determinant criterion** on the individually measured
  concentrations; Arrhenius reparameterised as k = X exp(-Y Ea), X = k0 exp(-Ea/(R T_av)),
  Y = (1/R)(1/T - 1/T_av), **T_av = 160 C**, so the 160 C column IS the estimated X. All five
  temperatures fitted simultaneously. **95 % HPD intervals.**
- **Reference temperature.** T_av = 160 C = **433.15 K** — identical to the registry's `T_REF_A_K`.
  No temperature transport is needed to compare these constants with the shipped ones.

## 3. Tables re-typed

### Table 1. "Model discrimination results for the proposed model based on Scheme 4 with (A) or without (B) the reversible isomerisation reaction from glucose to fructose (k6)."

| Model | p | SS | n | AICc | ΔAICc | PPB | PPS |
|---|---:|---|---:|---:|---:|---:|---:|
| A | 18 | 1.14 x 10^5 | 700 | 18182.7 | 48.2 | 74.5 | 0.05 |
| B | 16 | 1.13 x 10^5 | 700 | 18134.6 | 0.0 | 73.2 | 0.95 |

Footnote: "p (number of parameters); SS (residual sum of squares); n (number of data points
including the replicates); AICc (corrected Akaike criterion); ΔAICc (AICc difference taking the
smallest value as reference); PPB (Log10 of posterior probability); and PPS (normalised posterior
probability share)." Verdict: "The results from our test, shown in Table 1, support the model of
Scheme 4B without the reversible sugar isomerisation reaction." **Table 2 therefore reports model
4B.** Note on signs: PPB is defined as a log10 posterior probability and must be negative; the text
layer of this PDF drops leading minus signs in places (the same drop is visible in the
De Vleeschouwer Part II table, where "EaC 6.7" is printed "-6.7" in the neighbouring column), so
read PPB as **-74.5 and -73.2**. AICc and ΔAICc are internally consistent as printed
(18182.7 - 18134.6 = 48.1 ~ 48.2) and are positive.

### Table 2. "Estimates of rate constants (k) at 120, 140, 160, 180 and 200 C and activation energies (Ea) ±95 % highest posterior density (HPD) interval as found by kinetic modelling for the proposed kinetic model presented in Scheme 4B."

**Unit footnotes exactly as printed: `* (10^-3 min^-1)` and `** (10^-3 l mmol^-1 min^-1)`.**
Only k2 carries the double asterisk.

| k | 120 C | 140 C | 160 C | 180 C | 200 C | Ea (kJ/mol) |
|---|---|---|---|---|---|---|
| k1 * | 0.55 ± 0.1 | 1.4 ± 0.2 | 3.1 ± 0.4 | 6.5 ± 1 | 13 ± 3 | 61 ± 8 |
| **k2 \*\*** | **0.020 ± 0.01** | **0.078 ± 0.03** | **0.27 ± 0.06** | **0.86 ± 0.2** | **2.4 ± 0.6** | **93 ± 12** |
| k3 * | 1.9 ± 2 | 3.1 ± 3 | 5.0 ± 3 | 7.8 ± 3 | 12 ± 4 | 35 ± 27 |
| k5 * | 30 ± 37 | 55 ± 48 | 95 ± 52 | 157 ± 51 | 246 ± 82 | 40 ± 26 |
| k7 * | 1.7 ± 0.7 | 5.2 ± 1 | 14 ± 3 | 36 ± 7 | 84 ± 21 | 75 ± 10 |
| k8 * | 3.9 ± 2 | 8.3 ± 3 | 16 ± 5 | 30 ± 10 | 53 ± 23 | 50 ± 16 |
| k9 * | 2.0 ± 1 | 7.1 ± 2 | 22 ± 5 | 62 ± 13 | 159 ± 48 | 84 ± 14 |

**k4 and k6 are absent from the table because they were removed from the model**: k4 = acrylamide
degradation, deleted because its Ea came out negative; k6 = glucose -> fructose back-isomerisation,
deleted by the Table 1 discrimination.

**Which step is which.** Schemes 3 and 4 are images. The running text names five of the nine
constants directly and unambiguously:

- **k2** — "the formation of the Schiff base from fructose and asparagine (k2)": **Fru + Asn ->
  Schiff base**, and the ** footnote makes it second order.
- **k3** — "the formation of acrylamide from the Schiff base (k3)": **Schiff base -> acrylamide**.
- **k4** — "the activation energy for the breakdown of acrylamide was negative ... the degradation
  route of acrylamide was omitted": **acrylamide -> X1, DELETED**.
- **k6** — "removing the reversible isomerisation reaction of glucose to fructose (k6)":
  **Glc -> Fru, DELETED**.
- **k7 (in the earlier Scheme 1 numbering)** — "The model of reaction network presented in Scheme 1
  resulted in negative values ... for the formation of the Schiff base from glucose and asparagine
  (k6) and for the isomerisation reaction of glucose to fructose (k7)." **Careful: the numbering
  changed between Scheme 1 and Scheme 4.** In Scheme 1, k6 = Glc + Asn -> Schiff and k7 = Glc ->
  Fru; in Scheme 4, k6 = Glc -> Fru. The Scheme-1 step "Glc + Asn -> Schiff base" was dropped
  entirely ("the formation of the Schiff base at pH 5.5 from glucose and asparagine is not
  significant in a fructose-asparagine reaction system").

That leaves **k1, k5, k7, k8, k9** in the Scheme-4B numbering. My reading, argued:

- **k1 = fructose -> glucose (isomerisation)**, by elimination: it is the only isomerisation left
  after k6 was deleted, and the text says "The isomerisation of fructose to glucose was three times
  lower than was the isomerisation of glucose to fructose at pH 6.8 (Knol et al., 2005)" —
  3.1e-3 here against Knol 2005's k2 = 50.1e-3 at 160 C is 16x, not 3x, so the "three times"
  is either at a different temperature or loose; the identity is nonetheless forced.
- **k5 = Schiff base -> melanoidins**, and this one is confirmed by an independent coincidence:
  Ea(k5) = **40 ± 26** here against Knol 2005's k5 Ea = **40.1 ± 5.0**, the same step in the same
  slot of the same laboratory's earlier network, agreeing to one decimal place. Rates 95e-3 here
  vs 48.7e-3 there at 160 C, a factor 2.
- **k7 = fructose -> acetic acid** and **k9 = glucose -> formic acid**. The text fixes the pairing
  of substrate to product — "the degradations of fructose and glucose were included by the pathways
  leading to acetic acid and formic acid, respectively" — but not which k is which. **The magnitudes
  decide it (mine):** fructose is at 100 mmol/L and glucose only ever appears as a minor
  isomerisation product, yet acetic acid reaches "a maximum yield of 120 % (% mmol/mmol fructose)"
  while formic acid "never exceeded more than 15 % of the degraded D-fructose". A first-order
  constant acting on the small glucose pool must be the larger of the two to produce any formic
  acid at all; k9 = 22e-3 > k7 = 14e-3 at 160 C, so **k9 is the glucose (formic) route and k7 the
  fructose (acetic) route.** This matches the reading already carried in
  `k3_final_parameter_inventory.md` §A.1 ("acetic acid formation 75 ± 10; formic acid formation
  84 ± 14").
- **k8 = asparagine -> X2 (unidentified breakdown products)**, by elimination — "The degradation of
  asparagine was incorporated by including the putative formation of unknown reaction products
  (X2)."

**Read Schemes 3 and 4B from the page images before k1, k5, k7, k8 or k9 enters a registry.** k2
and k3 do not depend on any of this.

### Numbers printed in the running text (everything else is figure-only)

| quantity | value | where |
|---|---|---|
| **acrylamide yield** | "reaching almost **3 % of the initial asparagine concentration**" (against 1.5 % for glucose at pH 6.8 in Knol 2005) | Results 3.1; Conclusions |
| complete loss of fructose | "almost a complete loss" at 32 min (180 C) and 16 min (200 C) | Results 3.1, Fig. 2A |
| complete loss of asparagine | at 32 min (180 C) and 8 min (200 C) | Results 3.1, Fig. 2B |
| mannose | not formed; glucose formed | Results 3.1, Fig. 2C |
| acrylamide time shape | rises then "a slow decrease or steady state" at 160, 180, 200 C; the maximum coincides with asparagine (and partly fructose) reaching zero | Results 3.1, Fig. 2D |
| formic acid yield | "never exceeded more than 15 % of the degraded D-fructose"; three times as much as reported for glucose-glycine at 100 C / 4 h | Results 3.1 |
| **acetic acid yield** | "reaching a maximum yield of **120 % (% mmol/mmol fructose)**" | Results 3.1 |
| acetic > formic | "acetic acid was always formed in higher concentrations than was formic acid, regardless of the temperature" | Results 3.1 |
| pH | buffer 5.5 -> **5.3 on dissolving the reactants**; falls further on heating, maximum decrease at 32 min (180 C) and 16 min (200 C); at 200 C the final pH exceeds the initial pH; measured after cooling | Results 3.1, Fig. 2G |
| mass balance (% of initial fructose) | ~80 % at 120 C; ~100 % at 140 C rising to ~120 % at 64 min; ~100 % in the first 4-8 min at 160-200 C, then above 100 %, falling to ~85 % (180 C, 64 min) and ~65 % (200 C, 16-64 min) | Results 3.1, Fig. 3 |
| eps (melanoidin) | 282 L mol^-1 cm^-1, a **glucose/asparagine** value | Methods 2.7 (from Leong 1999) |
| stoichiometric rule imposed | "1 mol of fructose or glucose can, theoretically, be broken down to **3 mol** of acetic or formic acid" | Results 3.3 |
| k2 vs Knol 2005 | "lower than that at pH 6.8" (0.27e-3 here vs 1.45e-3 there at 160 C = **5.4x lower**, mine) | Results 3.3 |
| k3 vs Knol 2005 | "at 120 C, about a factor of 2 higher than that at pH 6.8" — **this does not reproduce**: 1.9e-3 here against Knol 2005's k4 = 0.176e-3 is 10.8x (mine). See Flags 3 | Results 3.3 |
| Ea(k3) vs Knol 2005 | "the activation energy for the formation of acrylamide from the Schiff base was almost a factor of 3 lower" — **this DOES reproduce**: 94.4 (Knol 2005 k4) / 35 (here) = 2.7 (mine), and it independently confirms that Knol 2005's k4 is the acrylamide-forming step | Results 3.3 |
| highest barrier in the network | "the Ea for the Schiff base formation is the highest in this network" — 93 ± 12, confirmed against the table | Results 3.3 |

**Concentration-time data: FIGURE-ONLY.** Figs. 2A-H and 5A-G carry every concentration, every
error bar, and the whole pH trajectory. Per house rule none is typed as a number.

### Arithmetic on the printed constants (all mine)

**1. The second-order unit, converted once.** 1 L mmol^-1 = 1000 L mol^-1, so a value of
x in units of 10^-3 L mmol^-1 min^-1 is **numerically x M^-1 min^-1** and **x x 10^-3
L/(mmol*min)** in the registry's own unit. k2 at 160 C: 0.27e-3 L mmol^-1 min^-1 =
**0.27 M^-1 min^-1 = 2.7e-4 L/(mmol*min)**.

**2. Two-point Arrhenius cross-check.** Ea = R ln(k200/k120)/(1/393.15 - 1/473.15) with the
denominator 4.302e-4: k1 gives 8.314e-3 x ln(13/0.55)/4.302e-4 = **61.1** (printed 61); k2
**92.4** (printed 93); k3 **35.5** (printed 35); k5 **40.7** (printed 40); k7 **75.2**
(printed 75); k8 **50.4** (printed 50); k9 **84.6** (printed 84). **Every row reproduces to
better than 1 %**, so each row carries two free parameters presented as six numbers, exactly as in
Knol 2005.

**3. Ratios inside the table at 160 C.** k5/k3 = 95/5.0 = **19.0** — the Schiff base goes to
browning 19 times more often than to acrylamide. Knol 2005's corresponding ratio (k5/k4) is
**19.2**. The two Wageningen pots, different sugar, different pH, different concentration, give
the same partition to within 1 %. That is the strongest internal result in the two papers and it is
directly relevant to `k_int1_mel`, the partition the shipped lane FITS because De Vleeschouwer's
k_M is refused. k9/k7 = 22/14 = **1.57**.

**4. What the 3 % yield is in absolute terms (mine).** 3 % of 100 mmol/L asparagine = 3.0 mmol/L
acrylamide = 3.0e-3 x 71.08 g/mol = **0.213 g/L = 2.1 x 10^5 ug/kg (ppb)** at unit density. Knol
2005's 1.5 % of 200 mmol/L is the same 3.0 mmol/L and the same 2.1 x 10^5 ppb. This is the number
`k3_final_parameter_inventory.md` §A.2 uses as the arbiter of the shipped lane's 480x two-lane
contradiction; it is confirmed here from the primary text.

**5. The elimination that is not here.** With k3 = 5.0e-3 min^-1 at 160 C feeding acrylamide and
**no sink at all**, this model's acrylamide can only rise or plateau. Fig. 2D shows "a slow
decrease or steady state" at 160-200 C and the fit "underestimated" acrylamide from 4 to 32 min at
those temperatures. If Knol 2005's own k6 (88.1e-3 min^-1 at 160 C, Ea 85.1) had been carried into
this pot unchanged, the acrylamide half-life would be 7.9 min and the 64-minute runs would show a
steep fall, which they do not. **Either the elimination is genuinely much slower at pH 5.5 with
fructose than at pH 6.8 with glucose, or one of the two fits is absorbing the difference
elsewhere.** The paper does not resolve this and neither can this dossier; it is recorded as the
sharpest open question the three Knol papers raise about `k_acr_dp`.

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`).** `acrylamide` is keyed. **Fructose, glucose,
asparagine, the Schiff base, melanoidins, acetic acid, formic acid, X1 and X2 are NOT in the
registry** (checked: 75 ids; the only acid present is `chlorogenic_acid`). Shared conditions for
every row: **100 mmol/L fructose + 100 mmol/L asparagine, 0.1 M phosphate, INITIAL pH 5.5 falling
to 5.3 on dissolution and drifting further on heating, dilute aqueous a_w ~ 1.0, sealed stainless
steel tube, heating block, 120-200 C, 0-64 min, duplicate/triplicate, heat-up NOT corrected,
T_ref = 160 C = 433.15 K.**

| step | quantity | value | unit as printed | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| **Fru + Asn -> Schiff base** | **k2 at 120/140/160/180/200 C** | 0.020 ± 0.01 / 0.078 ± 0.03 / **0.27 ± 0.06** / 0.86 ± 0.2 / 2.4 ± 0.6 | **10^-3 L mmol^-1 min^-1** (footnote \*\*) = M^-1 min^-1 numerically | aqueous, initial pH 5.5 | **second order**, mass action in [Asn][Fru] | Table 2 p. 1056 | measured_rate — **the corpus's only second-order Maillard-initiation constant from outside Leuven** |
| " | Ea(k2) | **93 ± 12** | kJ/mol | 120-200 C | — | Table 2 | measured_barrier (**the highest in this network**) |
| **Schiff base -> acrylamide** | **k3** | 1.9 ± 2 / 3.1 ± 3 / **5.0 ± 3** / 7.8 ± 3 / 12 ± 4 | 10^-3 min^-1 | as above | first order in the unmeasured intermediate | Table 2 | measured_rate — **HPD >= estimate at 120 and 140 C** |
| " | Ea(k3) | 35 ± 27 | kJ/mol | 120-200 C | — | Table 2 | measured_barrier — interval nearly spans zero; treat as weak |
| Fru -> Glc (isomerisation) | k1 | 0.55 ± 0.1 / 1.4 ± 0.2 / **3.1 ± 0.4** / 6.5 ± 1 / 13 ± 3 | 10^-3 min^-1 | as above | first order in fructose | Table 2 | measured_rate (step identity by elimination; Flags 2) |
| " | Ea(k1) | **61 ± 8** | kJ/mol | 120-200 C | — | Table 2 | measured_barrier — a trunk number, `k3_final_parameter_inventory.md` recommends FITTING it |
| Schiff base -> melanoidins | k5 | 30 ± 37 / 55 ± 48 / **95 ± 52** / 157 ± 51 / 246 ± 82 | 10^-3 min^-1 | as above | first order | Table 2 | measured_rate — **HPD > estimate at 120 and 140 C**; identity corroborated by Ea agreeing with Knol 2005's k5 to 40 vs 40.1 |
| " | Ea(k5) | 40 ± 26 | kJ/mol | 120-200 C | — | Table 2 | measured_barrier, weak |
| Fru -> acetic acid (x3 stoichiometry imposed) | k7 | 1.7 ± 0.7 / 5.2 ± 1 / **14 ± 3** / 36 ± 7 / 84 ± 21 | 10^-3 min^-1 | as above | first order in fructose | Table 2 | measured_rate (identity argued in section 3; Flags 2) |
| " | Ea(k7) | **75 ± 10** | kJ/mol | 120-200 C | — | Table 2 | measured_barrier — a trunk number |
| Asn -> X2 (unidentified) | k8 | 3.9 ± 2 / 8.3 ± 3 / **16 ± 5** / 30 ± 10 / 53 ± 23 | 10^-3 min^-1 | as above | first order in asparagine | Table 2 | measured_rate — **product not identified**; identity by elimination |
| " | Ea(k8) | 50 ± 16 | kJ/mol | 120-200 C | — | Table 2 | measured_barrier |
| Glc -> formic acid (x3 stoichiometry imposed) | k9 | 2.0 ± 1 / 7.1 ± 2 / **22 ± 5** / 62 ± 13 / 159 ± 48 | 10^-3 min^-1 | as above | first order in glucose | Table 2 | measured_rate (identity argued in section 3) |
| " | Ea(k9) | **84 ± 14** | kJ/mol | 120-200 C | — | Table 2 | measured_barrier — a trunk number |
| **acrylamide -> X1 (ELIMINATION)** | — | **NO VALUE. The step was fitted in Scheme 3, its activation energy came out NEGATIVE, and it was deleted from Scheme 4.** | — | — | — | Results 3.3 | **null result — high value as a null**, see section 1 |
| Glc + Asn -> Schiff base | — | **NO VALUE. Fitted negative in Scheme 1 and removed** ("not significant in a fructose-asparagine reaction system") | — | — | — | Results 3.3 | null result |
| Glc -> Fru (back-isomerisation) | k6 | **removed** by model discrimination (ΔAICc 48.2; PPS 0.95 vs 0.05) | — | — | — | Table 1 | within_study_ratio (structural) |
| browning-vs-acrylamide partition | k5/k3 at 160 C | **19.0** (Knol 2005's k5/k4 = 19.2) | — | as above | — | derived (mine) | within_study_ratio — rests on the k5 identity |
| acrylamide yield | ~3 % of initial asparagine = 3.0 mmol/L = **2.1 x 10^5 ppb** | % / mmol/L / ppb | 0.1 M Asn, pH 5.5, fructose; the maximum over the run | — | Results 3.1, Conclusions (level from text; time and T from Fig. 2D) | level_only |
| acetic acid yield | up to **120 %** (mmol per mmol fructose) | % | 160-200 C maxima | — | Results 3.1 | level_only |
| formic acid yield | **<= 15 %** of the degraded fructose | % | all T | — | Results 3.1 | level_only |
| pH trajectory, mass balances, all concentration-time courses | — | — | — | — | Figs. 2A-H, 3A-E, 5A-G | **figure_only** |

### Can these be put on the same basis as the trunk's constants? Step by step.

**T_ref is 160 C on both sides and k2's unit is second order and footnoted**, so the comparison is
arithmetic, not interpretive. This is the only paper in the five where that is true without an
inference.

**(a) The condensation — a real cross-laboratory check, and it passes.**

| source | step | k at 160 C | in the registry's unit | Ea (kJ/mol) | matrix |
|---|---|---|---|---|---|
| De Vleeschouwer 2009 Part I, glucose column | Glc + Asn -> Int1 (`k_asn_glc`) | 1.70 ± 1.05 M^-1 min^-1 | 1.70e-3 L/(mmol*min) | 117.5 ± 25.2 | a_w 0.92 powder, ~2.7 mol/kg |
| De Vleeschouwer 2009 Part I, fructose column | Fru + Asn -> Int1 | 0.22 ± 0.38 M^-1 min^-1 | 2.2e-4 L/(mmol*min) | 149.1 ± 87.7 | a_w 0.92 powder |
| **Knol 2010 (this paper), k2** | **Fru + Asn -> Schiff** | **0.27 ± 0.06 M^-1 min^-1** | **2.7e-4 L/(mmol*min)** | **93 ± 12** | **aqueous 0.1 mol/L, initial pH 5.5** |
| Knol 2005, k3 (unit inferred) | Fru + Asn -> Schiff | 1.45 ± 0.44 M^-1 min^-1 | 1.45e-3 L/(mmol*min) | 102 ± 14 | aqueous 0.2 mol/L, pH 6.8 |
| Knol 2005, k1 (unit inferred) | Glc + Asn -> Schiff | 0.668 ± 0.13 M^-1 min^-1 | 6.68e-4 L/(mmol*min) | 57.6 ± 8.0 | aqueous 0.2 mol/L, pH 6.8 |

**Fructose + asparagine: 0.27 (Wageningen, water, pH 5.5) against 0.22 ± 0.38 (Leuven, powder,
a_w 0.92) — a factor of 1.23, well inside the Leuven interval (mine).** Barriers 93 ± 12 against
149.1 ± 87.7 — overlapping, though Leuven's is nearly unconstrained.
**Glucose + asparagine: 0.668 (Wageningen 2005, unit inferred) against 1.70 ± 1.05 (Leuven, the
shipped `k_asn_glc`) — a factor of 2.5, inside the Leuven interval (lower bound 0.65) (mine).**
Barriers 57.6 ± 8.0 against 117.5 ± 25.2 — **not overlapping** (65.6 against 92.3).
The pattern across every comparable step in this five-paper set is the same: **the rates agree at
T_ref, the barriers do not.**

**(b) The acrylamide-forming step.** k3 = 5.0e-3 min^-1 at 160 C against the shipped `k_int1_acr`
= 3.57e-3 min^-1 (De Vleeschouwer Part I glucose) and Knol 2005's k4 = 2.53e-3 — **a 2.0x spread
across two laboratories, two sugars and two pH values.** Barriers 35 ± 27 / 159.2 ± 29.5 /
94.4 ± 11 — a 4.5x spread with no two intervals overlapping except by way of k3's near-zero lower
bound. The **rate** at 160 C is transportable; **the barrier is the disputed quantity in this
lane, on every step, and no single source can settle it.**

**(c) The partition.** k5/k3 = 19.0 here and k5/k4 = 19.2 in Knol 2005. If a re-calibration wave
wants a prior for `k_int1_mel`, this is it: **the browning fate of the intermediate is ~19x the
acrylamide fate at 160 C, reproduced in two independent Wageningen pots.** It is a within-study
ratio in each paper and a cross-study agreement between them; it is not a Leuven number and it does
not depend on the a_w 0.92 matrix.

**(d) The trunk constants.** k1 (isomerisation, Ea 61 ± 8), k7 (acetic acid, Ea 75 ± 10) and k9
(formic acid, Ea 84 ± 14) belong to Module 4 rather than to the acrylamide lane. The inventory's
own recommendation is to **split this paper by step** — hold out the acrylamide steps (k2, k3) as
the cross-lab extrapolation test and fit the organic-acid and isomerisation steps — which this
dossier supports: k1/k7/k9 do not touch acrylamide and holding them out would orphan three
"missing chemistry" lanes.

**(e) What cannot be transported.** No water-activity variation. No competing amino acid. No real
food. **The pH is not a controlled variable but a drifting one** (5.5 -> 5.3 -> lower -> partly
back), so this paper cannot be paired with Knol 2005 to make a pH ratio: sugar, concentration and
pH all changed at once. Every concentration is figure-only, so **no benchmark row**, only
constants and the 3 % / 2.1e5 ppb yield anchor.

## 5. Flags

1. **There is no acrylamide elimination constant in this paper, by the authors' decision.** The
   step existed in Scheme 3, produced a negative activation energy, and was removed. Anyone
   citing "Knol 2010" for an acrylamide degradation rate is citing a step that does not exist in
   the published model. `parameters_acrylamide.py`'s refusal entry is confirmed.
2. **k1, k5, k7, k8 and k9 are assigned to steps by argument, not by the paper.** Schemes 3 and 4B
   are images and the text names only k2, k3, k4 and k6. Section 3 gives the arguments (elimination
   for k1 and k8; an Ea coincidence with Knol 2005 for k5; a magnitude-versus-substrate-pool
   argument for k7 and k9). **Read Schemes 3 and 4B from the page images before any of these
   enters a registry.** k2 and k3 are safe.
3. **The paper's own comparison of k3 against Knol 2005 does not reproduce.** "The rate constant
   for the formation of acrylamide from the Schiff base (k3) at pH 5.5 was, at 120 C, about a
   factor of 2 higher than that at pH 6.8 (Knol et al., 2005)": the printed values are 1.9e-3 here
   and 0.176e-3 there, a factor of **10.8** (mine). Note that k3's HPD at 120 C is ±2 on 1.9, so
   the statement may be a loose reading of an interval that nearly reaches zero. The *Ea*
   comparison in the same paragraph does reproduce exactly (94.4/35 = 2.7 ~ "a factor of 3"),
   which is what licenses the k4 identification in `knol2005_extraction.md`.
4. **The paper's "three times lower" isomerisation comparison also does not reproduce**: k1 =
   3.1e-3 at 160 C here against Knol 2005's k2 = 50.1e-3, a factor of 16 (mine). Sugar, pH and
   concentration all differ, so this is not necessarily an error, but the quoted factor is not the
   ratio of the printed 160 C constants.
5. **Four of the seven rows have an HPD at or above the estimate somewhere in the series**: k3 at
   120 C (1.9 ± 2) and 140 C (3.1 ± 3); k5 at 120 C (30 ± 37) and 140 C (55 ± 48); k8 at 120 C
   (3.9 ± 2). The authors say so: "The precision of some of the parameters is low, which indicates
   that the number of data points should increase to improve precision." **k2, k7 and k9 are the
   well-determined rows.**
6. **The heat-up is absorbed into the constants.** Fig. 1 shows the profile; the fit is isothermal
   from t = 0, and the earliest sample is at 1 min in a 10 mL steel tube. This biases the 180 and
   200 C constants downward and flattens the fitted barriers — the same defect as Knol 2005, and a
   candidate explanation for why the Wageningen barriers are systematically below the Leuven ones
   (Leuven integrates a 4 s thermocouple log).
7. **pH is uncontrolled, admitted, and hidden in the parameters.** Initial 5.5, immediately 5.3,
   then drifting by more than a unit and partly reversing at 200 C, measured only after cooling.
   `acrylamide_conditions.py`'s pH factors are indexed to initial pH for exactly this reason; if
   this paper is ever paired with Knol 2005 to test that factor, the test is confounded by the
   sugar and concentration changes as well.
8. **The 3 mol acid per mol sugar stoichiometry is imposed, not measured** ("we applied the rule
   that 1 mol of fructose or glucose can, theoretically, be broken down to 3 mol of acetic or
   formic acid"). k7 and k9 are conditional on it; a different stoichiometry rescales them.
9. **eps = 282 L mol^-1 cm^-1 is a glucose/asparagine value used on a fructose/asparagine pot**,
   and it is a room-temperature literature value the authors themselves say "was not explored at
   the high temperatures used in this study". k5 is fitted against that response and inherits the
   error, as does the k5/k3 partition. `k3_final_parameter_inventory.md` already flags eps as
   amino-acid-specific to a factor 2.3.
10. **Mass balances exceed 100 %** (to ~120 % at 140 C / 64 min) and fall to ~65 % at 200 C. The
    authors list the reasons (asparagine breakdown products counted against initial fructose; the
    melanoidin quantification). "Too high recoveries will undoubtedly lead to problems." Expect the
    constants to move on a refit with a corrected balance.
11. **Parameter count does not close, as in Knol 2005.** Table 2 prints 7 rate constants and 7
    activation energies = 14; Table 1 says model B has p = 16 and model A p = 18. The difference
    A - B = 2 is the removed k6 pair, consistent. The residual 2 parameters of model B are not
    identified anywhere. Same flag as `knol2005_extraction.md` Flags 4.
12. **Acrylamide by HPLC-UV at 210 nm, external standard, no internal standard, no MS**, in a pot
    that browns heavily and forms organic acids that also absorb at 210 nm (the organic-acid method
    quantifies at the same wavelength on a different column). This is the least selective
    acrylamide determination of the five papers alongside Knol 2005's.
13. **What this paper does NOT contain**: any acrylamide elimination constant; any glucose +
    asparagine condensation constant (fitted negative and removed); any water-activity or moisture
    variation; any competing amino acid; any real food; any identification of X1 or X2; any
    tabulated concentration; any supplementary material.
14. **What to request from the authors**: (i) the concentration-time data behind Figs. 2A-H (five
    temperatures x eight times x eight responses); (ii) the **pH trajectory of Fig. 2G as numbers**
    — it is the only measured pH-versus-time in this whole lane and would test
    `acrylamide_conditions.py`'s initial-pH construction directly; (iii) Schemes 3 and 4B at
    resolution, or the k-to-step mapping; (iv) the Scheme-3 acrylamide degradation estimate and its
    negative Ea, which is a *measurement* of unidentifiability and would be worth recording; (v)
    the identity of the two extra parameters in Table 1's p = 16.
15. **Registry gaps against `data/keys/compounds.yml`**: `acrylamide` present. **Fructose,
    glucose, asparagine, the Schiff base, melanoidins, acetic acid and formic acid absent.** The
    trunk's organic-acid lanes (which this paper's k7 and k9 are candidates for) have no keyed
    products at all.
