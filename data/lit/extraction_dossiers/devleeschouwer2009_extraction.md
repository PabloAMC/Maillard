# De Vleeschouwer 2009 Part II — EXTRACTION (equimolar asparagine + glucose + a THIRD amino acid — glutamine or cysteine — freeze-dried to a_w 0.92 at 4 C, sealed inox tubes in an oil bath at 120/140/160/180/200 C with a 4 s thermocouple log; six responses fitted simultaneously by the determinant criterion with the control system's constants held fixed)

### The paper the shipped cysteine channels come from, read from its own Table 3 for the first time — and its new content beyond the four siblings is threefold: it is the corpus's ONLY measured bimolecular acrylamide-scavenging constant (k_E2 = 49.36 ± 1.18 M^-1 min^-1, Ea 51.3 ± 1.5) and its ONLY measured competitor sugar-consumption constant; it REPRINTS Part I's whole control column as fixed values, which independently confirms the three trunk constants the registry ships and resolves Part I's "10^-3 mm^-1" unit typo; and its model discrimination positively RULES OUT the alternative mechanism, that cysteine works by inhibiting acrylamide formation.

**Source on disk:** `data/articles/devleeschouwer2009.pdf` (12 pp., **Food Chemistry 114 (2009)
535-546**).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/devleeschouwer2009.txt`, 1012 lines). **Table 3, the only parameter table in
the paper, was re-extracted with `pdftotext -f 8 -l 8 -raw` because `-layout` scatters its
superscripts onto neighbouring lines; the two readings agree cell for cell and the `-raw` reading
is what settles the units.** Table 1 came through clean. Table 2 came through clean. Schemes 1 and
2 and Figures 1-4 are images. **Every concentration-time datum in this paper is figure-only:
Figs. 1A-F (glutamine system) and 2A-F (cysteine system) carry acrylamide, glucose, asparagine, the
second amino acid, aspartic acid and melanoidins at five temperatures with parity-plot inserts, and
Fig. 3 (acrylamide yield vs time at 140 and 180 C) and Fig. 4 (melanoidins recalculated with the
glutamine extinction coefficient) are likewise images. There is not one absolute acrylamide
concentration printed anywhere in this paper.** No supplementary material.

**Repo status before this dossier.** This paper is the source of four of the six populated
`MEASURED_ACRYLAMIDE` rows in `src/kinetic_core/parameters_acrylamide.py` — `k_acr_cys`,
`k_cys_sink`, `k_cys_glc` and `k_asp_sink`, all cited as `_DV2_SOURCE`, "Table 3 p. 542, CYSTEINE
column" — and it is named in the module docstring, in `HOLDOUT_EXPOSURE_DISCLOSURE` and in the
`competition_mechanism` metadata string. It reached those rows through
`k3_final_parameter_inventory.md` sec. A.2 and `k1_kinetic_parameters.md` sec. 2c, and **had no
extraction dossier.** Its **GLUTAMINE column is a declared HOLD-OUT** (see Flags 1 before reading
section 3).

## 0. Identity

| field | value |
|---|---|
| Title | "Role of precursors on the kinetics of acrylamide formation and elimination under low moisture conditions using a multiresponse approach – **Part II: Competitive reactions**" |
| Authors | Kristel De Vleeschouwer, Iesel Van der Plancken, Ann Van Loey, Marc E. Hendrickx (corresponding) — Laboratory of Food Technology, LFORCE, Katholieke Universiteit Leuven, Kasteelpark Arenberg 22, B-3001 Heverlee, Belgium |
| Venue | **Food Chemistry 114 (2009) 535-546.** Received 9 April 2008, revised 21 August 2008, accepted 25 September 2008 |
| DOI | 10.1016/j.foodchem.2008.09.084 |
| Funding | IWT-Vlaanderen PhD grant (K. De Vleeschouwer); FWO postdoctoral fellowship (I. Van der Plancken) |
| Naming | subscripts from the printed Nomenclature: `F` formation, `E` elimination, `E2` the cysteine-dependent elimination, `INT` formation of Int1 from glucose and asparagine, `INT2` formation of melanoidins from glucose and the second amino acid, `M` Maillard (Int1 -> Int2), `B` browning, `C` caramelization, `Asp` aspartic acid formation, `X` consumption of aspartic acid through unidentified reactions, `Y` consumption of cysteine through unidentified reactions, `Glu` glutamic acid formation. `DP` = unidentified acrylamide degradation products |
| **Part I, also on disk** | `data/articles/devleeschouwer2009b.pdf` = **"...Part I: Effect of the type of sugar", Food Chemistry 114 (2009) 116-126** — a **different paper**, cited throughout Part II as "in press". See Flags 2; **no dossier is written for it here** |
| Companions on disk | `claeys2005_extraction.md` (the same laboratory's dilute-aqueous competitor panel), `knol2005_extraction.md`, `knol2009_extraction.md`, `knol2010_extraction.md` (the Wageningen second-laboratory set), `devleeschouwer2006_extraction.md` (pH 4/6/8), `devleeschouwer2007_extraction.md` and `devleeschouwer2008_extraction.md` (the water-activity series) |

## 1. Why it matters

The shipped acrylamide lane's competition mechanism is stated in
`acrylamide_registry_metadata`'s `competition_mechanism` string as

> "two named mass-action channels per competitor: consumption of the SHARED GLUCOSE pool (measured
> for cysteine as k_INT2), and Michael-acceptor scavenging of acrylamide (measured for cysteine as
> k_E2). There is no per-amino-acid yield multiplier in the registry and no place to put one."

**Both of those "measured" claims rest on one column of one table in this paper**, and until now
neither had been checked against the page. This dossier checks them and they hold: Table 3's
cysteine column prints k_E2ref = 49.36 ± 1.18 M^-1 min^-1 with Ea_E2 = 51.3 ± 1.5 kJ/mol, and
k_INT2ref = 0.26 ± 0.02 M^-1 min^-1 with Ea_INT2 = 30.3 ± 1.6, exactly as the registry carries them.
The other two rows, `k_cys_sink` (k_Yref = 0.35 ± 0.01, Ea_Y = 110.5 ± 8.5) and `k_asp_sink`
(k_Xref = 0.04 ± 0.01, Ea_X = 97.2 ± 8.3), likewise.

**What is genuinely new here, against the four siblings dossiered today.** Knol 2005 is an aqueous
glucose-asparagine network with no competitor at all; Knol 2010 is the same with fructose at
pH 5.5 and deletes its elimination step; Knol 2009 is potato crisps with no rate constant of any
kind; Claeys 2005 is the same Leuven laboratory's competitor panel but in a **dilute aqueous** pot
(0.01 mol/L, 0.05 M citrate, pH 6) where every competitor's effect is absorbed into a single lumped
first-order k_E per system. **None of them contains a bimolecular constant for the
acrylamide-scavenging step, and none contains a competitor's sugar-consumption constant.** This
paper is the only source in the corpus for either. Four further things it adds that no sibling has:

1. **A second printing of Part I's entire control column.** Table 3's cysteine column carries
   k_Fref = **3.57**, k_Eref = **0.10**, k_INTref = **1.70**, k_Mref = **1.23**, k_Bref = **3.90**,
   k_Aspref = **26.43**, Ea_INT = **117.5**, Ea_M = **105.7**, Ea_B = **180.3**, Ea_Asp = **105.4**,
   Ea_E = **113.2**, Ea_C = **-6.7**, all marked with footnote c ("Parameters estimated
   independently of the other parameters by fixing at their estimated values") — i.e. they are the
   control-system values from Part I, reprinted. **Three of them are what the registry ships as
   `k_asn_glc` (1.70 M^-1 min^-1, Ea 117.5), `k_int1_acr` (3.57e-3 min^-1, Ea 159.2) and
   `k_asn_asp` (26.43e-3 min^-1, Ea 105.4).** They can now be verified against a second printed
   table without opening the hold-out-adjacent Part I. They match.
2. **The unit typo in Part I is resolved by Part II.** `parameters_acrylamide.py` carries the flag
   `printed_unit_is_a_typo_10e-3_mm-1_means_min-1` on `k_int1_acr`, with the note "The source prints
   the unit as '(10^-3 mm^-1)', which the inventory identifies as a typo for min^-1". **Part II
   prints the same row header as `k_Fref (10^-3 min^-1)`, correctly, against the same value 3.57.**
   The inventory's inference is now a printed fact from the same authors, same journal, same volume,
   same table structure. The flag can be downgraded from an inference to a corroborated correction.
3. **A positive model-discrimination result that rules out the rival mechanism.** Section 3.3.2
   re-fits the cysteine system with the previously fixed Maillard, caramelization and
   acrylamide-formation parameters **freed**, obtaining a comparable fit and no significant change
   except a higher Ea_F, and concludes verbatim: "**The reduced acrylamide concentrations in the
   presence of cysteine can thus not be attributed to the suggested inhibiting effect of cysteine on
   the browning reactions and on acrylamide formation.**" The registry's decision to represent
   cysteine as *scavenging plus sugar competition* and not as a formation inhibitor is therefore not
   a modelling convenience — the paper tested the alternative and rejected it.
4. **The first evidence anywhere in the corpus that the SECOND ORDER of `k_acr_cys` is real.**
   The registry flags it `order_assumed_never_tested_by_the_source`, and within this paper alone that
   is correct — the order is imposed by Eqs. (15) and (19) on the strength of a qualitative
   observation in Stadler 2003. But putting this paper beside Claeys 2005 gives a **185-fold
   change in cysteine concentration in the same laboratory**, and the observed suppression scales
   with it: see section 3, arithmetic 4. A concentration-independent cysteine term would be wrong by
   a factor of ~660; the first-order-in-cysteine term is wrong by a factor of ~3.

The paper also carries the number that makes the whole cysteine lane worth having: **the acrylamide
yield in the cysteine system is 0.1-0.2 % of the control's** (Table 2), i.e. a >99 % reduction,
against Claeys' 54-69 % of control for the same amino acid in a dilute pot — the same effect,
270-fold larger, in the matrix that resembles food.

What this paper does NOT give: any absolute acrylamide concentration (all figure-only); any pH;
any water-activity variation (one a_w, 0.92); any initial reactant concentration (Flags 5); any
amino acid other than glutamine and cysteine; any identification of `DP`, `X` or `Y`; any
uncertainty on the fixed parameters.

## 2. Methods as they matter to a model

- **Composition.** "Model systems composed of three **equimolar** parts of reactants... The first
  and second part consisted of asparagine and glucose, whereas the third part contained either
  **L-glutamine or L-cysteine**." Reactants mixed in water, frozen with liquid nitrogen,
  freeze-dried (0.01 mbar, Alpha 2-4, Christ) to a powder. Chemicals: acrylamide 99.9 % (Bio-Rad),
  methacrylamide ≥99 % (Merck), butyramide ≥98 % (Fluka), L-asparagine / L-glutamine / L-cysteine /
  L-aspartic acid all ≥99.5 % (Fluka), D-glucose HPLC grade ≥99 % (Sigma or Fluka). **Hydration
  state of the asparagine and the glucose is not stated** (Flags 5).
- **Water activity, and the moisture that goes with it.** Powders divided into **1.15 g** portions
  in small open containers, dried over P2O5, then equilibrated in sealed jars over saturated
  **Sr(NO3)2 at 4 C for about 3 weeks** to **a_w 0.92** (Greenspan 1977). Total moisture by
  automated **Karl-Fischer titration** and printed in Table 1. The paper states the consequence
  itself: "Large differences can be observed between the moisture content of the model systems
  tested, which can be attributed to differences in hygroscopic properties of the different
  reactants used. **Consequently, the initial reactant concentration differs between the model
  systems tested.**"
- **Heating.** Hermetically sealed **custom-made inox reactor tubes, 8 mm x 100 mm**, filled in a
  closed dry-nitrogen environment (Captair pyramid, erlab) to prevent water resorption. Thermostated
  oil bath (UH2D, Grant) at **120, 140, 160, 180 and 200 C**. Sampling times "chosen dependent on
  the treatment temperature"; the figure time axes run to 130 min. Immediate ice-water cooling,
  10-fold dilution, storage at -40 C.
- **The thermal history is integrated, not assumed.** "During the heating and subsequent cooling
  phase of the samples, temperature was registered within the closed reactor tubes at regular time
  intervals (**4 s**) using thermocouples (type T)... The registered temperature-time profiles were
  inputted in the kinetic data analysis." This is the methodological difference from Knol 2005 and
  Knol 2010, both of which fit isothermally from t = 0.
- **Acrylamide.** GC-MS with **chemical ionisation**, method of De Vleeschouwer et al. 2006. Internal
  standards methacrylamide and butyramide are in the chemicals list. No detection limit is printed
  in this paper.
- **Amino acids.** EZ:faast kit (Phenomenex) by GC-MS, method of De Vleeschouwer et al. 2008a,
  quantified against the **norvaline** internal standard's m/z 158 ion. **"The concentration of
  cysteine was calculated based on the sum of the peak areas of the m/z 248 ion eluted at two
  different reaction times in order to account for both cysteine and cysteine spontaneously oxidised
  to cystine."** So `[Cys]` throughout is **cysteine + cystine**, and the registry's conditions
  string says so.
- **Sugars.** HPAEC-PAD (Dionex).
- **Browning.** Absorbance at 470 nm, converted to a melanoidin concentration by Lambert-Beer with
  **eps = 282 L/mol cm**, the coefficient for **asparagine-glucose** melanoidins, taken from Knol
  et al. 2005. The paper is explicit that this is wrong for a two-amino-acid system and says so
  three times; see Flags 6.
- **Fitting.** Differential equations per reaction step with the Arrhenius equation substituted in
  (Eq. 1, `k = k_ref exp(-(Ea/R)(1/T_ref - 1/T))`, R = 8.3145 J/mol K, **T_ref = 160 C**); numerical
  integration of each sample's registered temperature-time profile; **non-linear regression using
  the determinant criterion** (van Boekel 1996) in **Athena Visual Studio v11.0**; goodness of fit by
  scrutiny of residuals and parity plots. Uncertainties are **±95 % highest posterior density (HPD)
  intervals** (Table 3 footnote a).
- **The fixing strategy, which is the paper's whole logic.** "If the basic kinetic model proposed
  for the control system is consistent, then the different rate constants for the Maillard and
  caramelization reactions and the formation of aspartic acid should be independent of a change in
  the concentration of the reactants or reaction products due to addition of a second amino acid.
  Therefore, these reaction rate constants will be fixed to their value estimated for the control
  system." So `INT`, `M`, `B`, `C` and `Asp` are fixed from Part I in both columns; `X` is fixed in
  the glutamine column and freed in the cysteine column; `F` and `E` are freed in the glutamine
  column and fixed in the cysteine column. **Each printed number must be read together with its
  footnote letter, because roughly half of Table 3 is not a measurement made in this paper.**
- **Reference temperature.** T_ref = 160 C = 433.15 K, the registry's `T_REF_A_K`. No temperature
  transport is needed to place these constants beside the shipped ones.

### The two reaction networks, as printed in the differential equations

Scheme 1 (glutamine) and Scheme 2 (cysteine) are images, but the ODE systems are printed in full
and define the networks unambiguously.

**Glutamine, Eqs. (2)-(12), p. 541.** Asn + Glc -> Int1 (`k_INT`, second order); Int1 -> AA
(`k_F`); Int1 -> Int2 (`k_M`); Int2 + Asn -> browning (`k_B`, second order); Glc -> Int2
(caramelization, `k_C`); **Gln + Glc -> melanoidins in one lumped step (`k_INT2`, second order)**;
Gln -> Glu (`k_Glu`, deamidation, first order); Asn -> Asp (`k_Asp`); Asp -> X (`k_X`);
AA -> DP (`k_E`, first order). **Two printed typos in this ODE set are recorded in Flags 8.**

**Cysteine, Eqs. (13)-(23), p. 543.** The same, with Gln replaced by Cys and **two changes**:
`d[Cys]/dt = -k_INT2[Cys][Glc] - k_Y[Cys] - k_E2[Cys][AA]` and
`d[AA]/dt = k_F[Int1] - k_E[AA] - **k_E2[Cys][AA]**`. So the extra elimination is written as
**second order, first order in each of acrylamide and cysteine**, and the paper says why: "the
additional elimination reaction of acrylamide is **assumed** to be a second-order reaction,
depending on both acrylamide concentration and cysteine concentration, **based on the observation
of Stadler et al. (2003) that the removal of acrylamide in the presence of cysteine was enhanced
with increasing cysteine concentration.**" That sentence is the entire warrant for the order, inside
this paper, and it is the sentence behind the registry's
`order_assumed_never_tested_by_the_source` flag.

Note that `k_E2` is written to consume cysteine stoichiometrically (Eq. 15) as well as acrylamide
(Eq. 19) — a 1 : 1 Michael adduct, S-(2-carbamoylethyl)cysteine, which is the transformation string
the registry's `k_acr_cys` carries.

## 3. Tables re-typed

**Read Flags 1 first: the GLUTAMINE column of Tables 2 and 3 is a declared HOLD-OUT.**

### Table 1 (p. 536). "Initial moisture content of the equimolar asparagine-glucose-amino acid model systems equilibrated at a water activity of 0.92 (at 4 °C)."

Footnote a: "Standard error based on **two** independent measurements."

| Asn-Glc + Amino acid | Moisture content (%) |
|---|---|
| Control | **14.53 ± 0.08**^a |
| Glutamine | **9.58 ± 0.35**^a |
| Cysteine | **19.65 ± 0.12**^a |

### Table 2 (p. 539). "Effect of addition of an amino acid, other than asparagine, on the relative maximum acrylamide yield per mol initial asparagine concentration (%) in equimolar asparagine-glucose-amino acid model systems, equilibrated at an initial water activity of 0.92 (at 4 °C), heated at temperatures between 120 and 200 °C."

| T (°C) | Control | +Glutamine | +Cysteine |
|---|---|---|---|
| 120 | 100 | 267.2 | **0.2** |
| 140 | 100 | 180.1 | **0.2** |
| 160 | 100 | 132.4 | **0.2** |
| 180 | 100 | 132.4 | **0.1** |
| 200 | 100 | 120.0 | **0.1** |

**No absolute acrylamide concentration is printed anywhere in this paper.** Table 2 is entirely
relative and it is the paper's only tabulated response. The two 132.4 entries at 160 and 180 C are
identical to four figures (Flags 9).

### Table 3 (p. 542). "Estimated kinetic parameters based on multiresponse data describing acrylamide formation and elimination in an equimolar asparagine-glucose-amino acid model systems equilibrated at an initial water activity of 0.92 (at 4 °C), heated at temperatures between 120 and 200 °C."

Header: **T_ref = 160 °C**, columns "Asn-Glc-Amino acid": **Gln (Scheme 1)** and **Cys (Scheme 2)**.
Footnote a: "±95 % highest posterior density (HPD) interval." Footnote b: "Parameters not inserted
in the model." Footnote c: "Parameters estimated independently of the other parameters by fixing at
their estimated values." **There is no footnote `*`** — see Flags 7.

Units exactly as printed (confirmed in `-raw`; `-layout` drops the minus signs of the superscripts).

| parameter | unit as printed | **Gln** (Scheme 1) | **Cys** (Scheme 2) |
|---|---|---|---|
| k_Fref | (10^-3 min^-1) | 8.05 ± 0.90^a | **3.57**^c |
| k_Eref | (min^-1) | 0.36 ± 0.05 | **0.10**^c |
| k_E2ref | (M^-1 min^-1) | –^b | **49.36 ± 1.18** |
| k_INTref | (M^-1 min^-1) | 1.70^c | **1.70**^c |
| k_INT2ref | (M^-1 min^-1) | 0.01 ± 0.00 | **0.26 ± 0.02** |
| k_Mref | (min^-1) | 1.23^c | **1.23**^c |
| k_Bref | (M^-1 min^-1) | 3.90^c | **3.90**^c |
| k_Cref | (10^-3 min^-1) | Indeterminate^c | **Indeterminate**^c |
| k_Aspref | (10^-3 min^-1) | 26.43^c | **26.43**^c |
| k_Xref | (min^-1) | Indeterminate^c | **0.04 ± 0.01** |
| k_Gluref | (min^-1) | 0.62 ± 0.33 | –^b |
| k_Yref | (min^-1) | –^b | **0.35 ± 0.01** |
| Ea_F | (kJ/mol) | 124.1 ± 9.3 | **159.2*** |
| Ea_E | (kJ/mol) | 92.4 ± 12.0 | **113.2**^c |
| Ea_E2 | (kJ/mol) | –^b | **51.3 ± 1.5** |
| Ea_INT | (kJ/mol) | 117.5^c | **117.5**^c |
| Ea_INT2 | (kJ/mol) | 13.2 ± 4.3 | **30.3 ± 1.6** |
| Ea_M | (kJ/mol) | 105.7^c | **105.7**^c |
| Ea_B | (kJ/mol) | 180.3^c | **180.3**^c |
| Ea_C | (kJ/mol) | -6.7^c | **-6.7**^c |
| Ea_Asp | (kJ/mol) | 105.4^c | **105.4**^c |
| Ea_X | (kJ/mol) | 668.9^c | **97.2 ± 8.3** |
| Ea_Glu | (kJ/mol) | 35.9 ± 6.8 | –^b |
| Ea_Y | (kJ/mol) | –^b | **110.5 ± 8.5** |

The paper adds, in the text, that **k_INT2ref and Ea_INT2 "are thus only apparent values and
therefore they are indicated in grey in Table 3"**, along with the two `Glu` parameters. The grey
shading does not survive text extraction; the qualification is recorded here in its place, and the
registry already carries it as `apparent_only_authors_say_not_comparable_to_k_INT`.

### Numbers printed in the running text (everything else in this paper is figure-only)

| quantity | value | where |
|---|---|---|
| acrylamide level, cysteine system vs control | "the net acrylamide concentration of the system with cysteine is **more than 1000 times lower** than the concentrations measured for the control system and the system with additional glutamine" | Results 3.1, p. 537 |
| glutamine, effect on the yield | "the maximum is not only higher when glutamine is added..., but is also **attained earlier**"; the lag phase is "even shorter" at 120 and 140 C | Results 3.1, pp. 537-538 |
| glucose consumption | "the glucose concentration of only **a few of the samples exceeds 5 %** of the initial concentration due to its very rapid consumption" | Results 3.3.1, p. 542 |
| relative decrease rates | glucose faster than asparagine; asparagine slower than the added amino acid in both systems; **glutamine is consumed faster than cysteine** | Results 3.2, p. 541 |
| glutamic acid | detected, but "the glutamic acid concentrations measured were... so small that **no clear trend** as a function of time and temperature could be deduced" | Results 3.2, p. 541 |
| aspartic acid yield | "comparable for the control system and the system with cysteine added and is **slightly lower** for the system with glutamine added" (data not shown) | Results 3.2, p. 541 |
| **apparent melanoidin yield per mol glucose, relative to control** | glutamine system **111 %**; cysteine system **33 %** (data not shown) | Results 3.2, p. 541 |
| **eps for glutamine-glucose melanoidins** | **498 L/mol cm** (Leong 1999) against **282 L/mol cm** for asparagine-glucose | Results 3.2, p. 541 |
| eps for cysteine-glucose melanoidins | **none exists**; Ashoor & Zent (1984) give the colour intensity of cysteine-glucose browning as "one of the lowest" and "**about 25 %** of the intensity measured for an asparagine-glucose mixture" | Results 3.2, p. 541 |
| melanoidin recalculation | using the glutamine-glucose eps, "**none of the kinetic parameters estimated were changed significantly**... This proves that the kinetic parameters describing acrylamide formation and elimination reactions are independent of the extinction coefficient used" | Results 3.3.1, p. 542 |
| **the rival mechanism, tested and rejected** | freeing the fixed Maillard, caramelization and formation parameters in the cysteine system gave a comparable fit with no significant differences except a higher Ea_F: "**The reduced acrylamide concentrations in the presence of cysteine can thus not be attributed to the suggested inhibiting effect of cysteine on the browning reactions and on acrylamide formation.**" | Results 3.3.2, p. 544 |
| k_E2 vs k_E, in the authors' words | "The rate constant of the cysteine-dependent acrylamide elimination reaction (k_E2ref) is **quite high** as compared to rate constant of the basic acrylamide elimination reaction (k_Eref), **irrespective of its unit**. The temperature dependence of k_E2ref... is **about half** of the activation energy of the basic acrylamide elimination reaction" | Results 3.3.2, p. 544 |
| k_INT2, Cys vs Gln | the cysteine value "is significantly higher than the corresponding value estimated for the asparagine-glucose-glutamine system, but both values are nevertheless **in the same order of magnitude**"; the two Ea_INT2 are **not significantly different** from each other but **are** significantly lower than Ea_INT | Results 3.3.2, p. 544 |
| the cysteine sink `Y` is a lump | "the nature of the reaction consuming cysteine (indicated with the symbol 'Y') is **not specified** and this could even represent a **group of reactions** with each a different reaction rate and activation energy" | Results 3.3.2, p. 544 |
| browning fit, cysteine system | "The model **strongly overestimates** browning for all the temperature-time conditions applied" | Results 3.3.2, p. 544 |
| browning fit, glutamine system | melanoidins "generally **underestimated** by the model" | Results 3.3.1, p. 542 |
| why glutamine's k_F was freed | fixing it gave a fit that "was, however, **inadequate**" | Results 3.3.1, p. 542 |
| glutamine cannot itself make acrylamide | glutamine "can result in only **negligible** amounts of acrylamide in the presence of a sugar" (Leufvén & Lingnert 2003; Mottram 2002), and the purity ≥99.5 % excludes contamination as the explanation | Results 3.1, p. 538 |

**Figure-only in this paper (not typed as numbers).** Figs. 1A-F and 2A-F: acrylamide, glucose,
asparagine, the second amino acid, aspartic acid and melanoidins against time at five temperatures,
each with a parity-plot insert. **Every concentration in this paper, including every initial
concentration, is in those figures only.** Fig. 3: acrylamide yield per mol initial asparagine
against time for the control and glutamine systems at 140 and 180 C. Fig. 4: melanoidins
recalculated with the glutamine-glucose eps. Schemes 1 and 2: the two reaction networks.

### Arithmetic on the printed numbers (all mine)

**1. The registry's four rows reproduce exactly.** `k_acr_cys` = 49.36 M^-1 min^-1 converted by
`per_molar_to_per_mmol` is 4.936e-2 L mmol^-1 min^-1, with Ea 51.3 ± 1.5. Relative standard error
1.18/49.36 = **2.39 %**, which is the "2.4 % RSE" the registry calls the corpus's tightest
parameter — confirmed from the primary table. `k_cys_glc` = 0.26 M^-1 min^-1 -> 2.6e-4 L mmol^-1
min^-1, Ea 30.3 ± 1.6. `k_cys_sink` = 0.35 min^-1, Ea 110.5 ± 8.5. `k_asp_sink` = 0.04 min^-1,
Ea 97.2 ± 8.3. All four match the page.

**2. Initial concentrations, which the paper never prints, derived from Table 1's moisture.**
Taking the printed moisture as w/w of the equilibrated powder, the equimolar stoichiometry, a
density of 1 kg/L and anhydrous molecular weights (Asn 132.12, Glc 180.16, Gln 146.15, Cys 121.16):

| system | moisture | solids | Σ MW per equimolar set | **mol/kg of each reactant** |
|---|---|---|---|---|
| control (Asn + Glc) | 14.53 % | 854.7 g/kg | 312.28 | **2.74** |
| + glutamine | 9.58 % | 904.2 g/kg | 458.43 | **1.97** |
| + cysteine | 19.65 % | 803.5 g/kg | 433.44 | **1.85** |

With asparagine as the monohydrate instead the three become 2.59, 1.90 and 1.78 — so the **ratios**
are robust: **the control system is 1.36-1.39x more concentrated in asparagine and glucose than the
two competitor systems.** This is the quantitative content of the paper's own sentence "the initial
reactant concentration differs between the model systems tested", and it matters twice over
(Flags 4): the model does integrate each sample's measured initial concentrations, so the ODEs are
right, but **k_INT, k_M, k_B and k_Asp were fixed at values estimated at a 1.4x different
concentration**, and k_INT and k_B are the two second-order steps where a lumped apparent constant
is most likely to drift with concentration. The registry's `_DV2_CONDITIONS` string says "~3 mol/L";
**for the cysteine system the better number is ~1.85 mol/L**, and the difference is a factor of 1.6
on every bimolecular flux computed from it.

**3. The suppression the cysteine channel predicts, against the suppression Table 2 measures.**
At T_ref = 160 C, with [Cys]_0 = 1.854 mol/L from row 2 above, the cysteine elimination flux
coefficient is k_E2 [Cys] = 49.36 x 1.854 = **91.5 min^-1**, against the base k_E = **0.10 min^-1**
— a ratio of **915** at time zero. Table 2 measures a maximum-yield suppression of 1/0.002 = **500**
at 120-160 C and 1/0.001 = **1000** at 180-200 C. **Same order of magnitude, and it should be
slightly below the t = 0 estimate**, because cysteine is itself consumed fast: at 160 C its two
printed sinks give k_Y + k_INT2 [Glc] = 0.35 + 0.26 x 1.854 = 0.832 min^-1, a half-life of
**0.83 min**. So the paper's headline ">99 %" and its fitted constants are mutually consistent, which
is a real internal check that no sibling paper allows. The *temperature trend* cannot be checked:
Table 2's cysteine column carries **one significant figure**, and 0.2 -> 0.1 is not enough
resolution to test Ea_E2 = 51.3 against Ea_E = 113.2.

**4. The second-order form of `k_acr_cys`, tested across a 185-fold concentration change.**
This is the most useful derived result in the dossier, and it addresses the registry's own flag.
Claeys 2005 (same laboratory, same T_ref = 160 C, same GC-MS) fits **lumped** first-order
elimination constants: control k_Eref = 111.1e-3 min^-1, cysteine system k_Eref = 268.7e-3 min^-1,
in a **0.01 mol/L** aqueous pot. If the whole of that increase is the bimolecular cysteine channel,
then k_E2 = (268.7 - 111.1)e-3 / 0.01 = **15.8 M^-1 min^-1**. De Vleeschouwer II fits it directly at
**49.36 M^-1 min^-1** in a 1.85 mol/L powder. **A factor of 3.1 apart, across a 185-fold change in
cysteine concentration, a change of matrix from dilute aqueous at pH 6 to a freeze-dried powder at
a_w 0.92, and two entirely different fitting structures.** Note also that the two base elimination
constants agree almost exactly (111.1e-3 against 100e-3 min^-1), so the comparison is not resting on
a normalisation. **The alternative hypothesis — that cysteine's effect is concentration-independent
— predicts the same fractional suppression in both pots**, i.e. Claeys' 43 % reduction at 160 C
should also be what a_w 0.92 shows, against the measured >99 %: **wrong by a factor of about 660.**
I read this as substantial support for first order in cysteine, and I would put it at roughly 3 : 1
in favour rather than decisive, because the two pots differ in more than concentration and because
Claeys' Table 1 percentages are averages over nine sampling times while Table 2 here is a maximum
yield. **The registry's flag should be softened from "never tested" to "assumed by the source and
supported across two of its laboratory's own systems to within a factor of 3".**

**5. The two barriers cross over, and where.** Ea_E2 = 51.3 against Ea_E = 113.2 kJ/mol means the
cysteine channel's *relative* advantage falls as temperature rises. Using
k_E2[Cys]/k_E = 915 at 433.15 K and the two barriers, that ratio is **5.3e3 at 120 C** and
**2.1e2 at 200 C** (mine) — a 25-fold narrowing across the range, and it never reverses inside it.
Extrapolating the two Arrhenius lines to equality gives a crossover far above 200 C and outside any
licensed range; **the registry's statement that "the two channels CROSS OVER in temperature" is
directionally right but the crossing is not inside the measured window at this cysteine
concentration.** It would come inside the window at a much lower cysteine loading — which is the
same statement as row 4 above and is exactly the behaviour an order-2 term is for.

**6. Ea_C = -6.7 kJ/mol is negative, and the constant it belongs to is "Indeterminate".** A negative
activation energy on the caramelization step, attached to a rate constant the fit could not
determine in either column, is a fixed value inherited from Part I. **Neither k_C nor Ea_C is
usable**, and the registry correctly carries no caramelization row.

**7. Ea_X = 668.9 kJ/mol in the glutamine column is physically absurd** — roughly five times the
largest chemically meaningful Maillard barrier in the corpus, and attached to an "Indeterminate"
rate constant. The registry already refuses it and takes the cysteine column's 0.04 ± 0.01 /
97.2 ± 8.3 instead, with the note "The cysteine column's... is the one that was actually estimated".
**Confirmed from the page: 668.9 carries footnote c (fixed from Part I), 97.2 ± 8.3 carries an HPD
interval.**

## 4. Kinetic numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** `acrylamide` is keyed (CAS 79-06-1).
**Asparagine, glutamine, cysteine, cystine, aspartic acid, glutamic acid, glucose, melanoidins,
and the adduct S-(2-carbamoylethyl)cysteine are all absent** — the registry is a product/marker list
and carries no Maillard reactants (Flags 11). Every row below shares: equimolar Asn + Glc + a third
amino acid, freeze-dried powder equilibrated to **a_w 0.92 at 4 C**, ~1.85-1.97 mol/kg of each
reactant (mine, row 2 above), hermetically sealed inox tubes, **no pH** (a powder), 120-200 C,
4 s thermocouple log integrated into the fit, multiresponse determinant criterion,
**T_ref = 160 C = 433.15 K**.

| step | quantity | value | unit | conditions | reaction order as fitted | source location | evidence class |
|---|---|---|---|---|---|---|---|
| **AA + Cys -> S-(2-carbamoylethyl)cysteine** | **k_E2ref** | **49.36 ± 1.18** | **M^-1 min^-1** | cysteine system, a_w 0.92, 160 C | **second — first order in each of [AA] and [Cys], ASSUMED by the source (Eqs. 15, 19) on Stadler 2003's qualitative observation** | Table 3 p. 542, Cys column | **measured_rate** — the corpus's only bimolecular acrylamide-scavenging constant; 2.39 % RSE |
| " | **Ea_E2** | **51.3 ± 1.5** | kJ/mol | 120-200 C | — | Table 3 | **measured_barrier** — the lowest barrier in the table, less than half of Ea_E |
| **Cys + Glc -> melanoidins (LUMPED)** | k_INT2ref | **0.26 ± 0.02** | M^-1 min^-1 | as above | second | Table 3, Cys column | measured_rate — **carry the authors' words: "only apparent values"**, printed in grey, "not comparable" to k_INT because INT2 lumps initiation, intermediate and final steps while INT is initiation only |
| " | Ea_INT2 | **30.3 ± 1.6** | kJ/mol | " | — | Table 3 | measured_barrier (apparent, same caveat) |
| **Cys -> unidentified products** | k_Yref | **0.35 ± 0.01** | min^-1 | as above | first | Table 3, Cys column | measured_rate — **the product is not identified and the source says the step "could even represent a group of reactions"** |
| " | Ea_Y | **110.5 ± 8.5** | kJ/mol | " | — | Table 3 | measured_barrier |
| **Asp -> unidentified products** | k_Xref | **0.04 ± 0.01** | min^-1 | as above | first | Table 3, Cys column | measured_rate — **the only k_X in the corpus that was actually estimated**; both other columns print "Indeterminate" |
| " | Ea_X | **97.2 ± 8.3** | kJ/mol | " | — | Table 3 | measured_barrier |
| **Asn + Glc -> Int1** | k_INTref | **1.70** | M^-1 min^-1 | **fixed, footnote c** — the Part I control value | second | Table 3, both columns | **reprint of a Part I estimate** — no interval here; use only to *verify* the shipped `k_asn_glc` |
| **Int1 -> AA** | k_Fref | **3.57**, unit printed **(10^-3 min^-1)** | 10^-3 min^-1 | **fixed, footnote c** | first | Table 3, Cys column | **reprint** — and **the correct unit**, which is the new content (Flags 3) |
| **AA -> DP** | k_Eref | **0.10** | min^-1 | **fixed, footnote c** | first | Table 3, Cys column | reprint |
| **Int1 -> Int2** | k_Mref | **1.23** | min^-1 | fixed, footnote c | first | Table 3 | reprint |
| **Int2 + Asn -> browning** | k_Bref | **3.90** | M^-1 min^-1 | fixed, footnote c | second | Table 3 | reprint |
| **Asn -> Asp + NH3** | k_Aspref | **26.43** | 10^-3 min^-1 | fixed, footnote c | first | Table 3 | reprint |
| **Glc -> Int2 (caramelization)** | k_Cref | **Indeterminate** | — | fixed, footnote c | — | Table 3 | **not usable** |
| barriers reprinted from Part I | Ea_INT / Ea_F / Ea_E / Ea_M / Ea_B / Ea_C / Ea_Asp | **117.5 / 159.2 / 113.2 / 105.7 / 180.3 / -6.7 / 105.4** | kJ/mol | fixed, footnote c (Ea_F carries an **undefined** `*`) | — | Table 3, Cys column | **reprint** — verification of the shipped values; Ea_C is negative and unusable |
| **relative maximum acrylamide yield, cysteine** | **0.2 / 0.2 / 0.2 / 0.1 / 0.1** | % of control at 120/140/160/180/200 C | % | a_w 0.92 | — | Table 2 p. 539 | **within_study_ratio** — one significant figure; the paper's only tabulated response |
| moisture at a_w 0.92 | 14.53 ± 0.08 / 9.58 ± 0.35 / 19.65 ± 0.12 | % (control / Gln / Cys) | Karl-Fischer, n = 2 | — | Table 1 p. 536 | measured level |
| initial reactant concentrations | 2.74 / 1.97 / 1.85 | mol/kg | derived from Table 1 (mine) | — | — | **derived_assumption** — assumes w/w moisture, density 1, anhydrous MWs |
| eps, melanoidins from Gln + Glc | 498 | L/mol cm | — | — | Results p. 541, from Leong 1999 | **level_only** — borrowed, not measured here |
| eps, melanoidins from Asn + Glc | 282 | L/mol cm | — | — | Methods p. 537, from Knol 2005 | **level_only** — borrowed |
| cysteine-glucose browning intensity | ~25 % of asparagine-glucose | — | Ashoor & Zent's conditions, not these | — | Results p. 541 | **level_only** — borrowed, and from a different experiment |
| apparent melanoidin yield per mol glucose vs control | 111 % (Gln), 33 % (Cys) | % | a_w 0.92 | — | Results p. 541, "data not shown" | within_study_ratio (weak — the underlying values are not shown and the eps is wrong for both systems) |
| every concentration-time course, every parity plot, the acrylamide yield curves | — | mM vs min | — | — | Figs. 1A-F, 2A-F, 3, 4 | **figure_only** |

### Can these be put beside the shipped constants? Yes, and here is the arithmetic that says so.

**T_ref matches exactly** (160 C in both), so no temperature transport is needed. **a_w matches
exactly** for the four cysteine rows, which is why the registry marks them
`licensed_at_measurement_aw_only` — with the one exception of `k_acr_cys`, whose transfer licence is
widened to `licensed_to_the_thiol_michael_acceptor_family`. That widening is not licensed *by this
paper*, which measures one thiol at one water activity; it comes from the inventory's verdict, and
this dossier does not add to it.

**What the cysteine column proves about the trunk.** The whole point of fixing `INT`, `M`, `B`, `C`
and `Asp` is a consistency test, and the paper states the conclusion: "it can be assumed that the
basic model, used to describe the reactions in an asparagine-glucose model systems under the same
reaction conditions, remains unaffected by addition of an equimolar concentration of cysteine and is
thus consistent." **So the shipped `k_asn_glc`, `k_int1_acr` and `k_asn_asp` have survived a
consistency test in a third system**, which no sibling paper provides. The test is weaker than it
looks — see Flags 4 — because the fixed constants were estimated at a 1.4x higher reactant
concentration and because the browning response, the one the fixed `k_B` and `k_C` act on, is
badly fitted in both columns.

**What the paper does NOT let the repository do.** It cannot supply a pH term (a powder has no pH),
a water-activity exponent (one a_w), any competitor other than cysteine on the FIT side, any
benchmark row (no absolute concentration is printed), or any identification of `DP`, `X` or `Y`.
And it cannot, on its own, test the order of the scavenging step — that took a second paper
(section 3, arithmetic 4).

## 5. Flags

1. **The GLUTAMINE column of Tables 2 and 3 is a DECLARED HOLD-OUT and this dossier prints it in
   full.** `HOLDOUT_EXPOSURE_DISCLOSURE` in `parameters_acrylamide.py` already records that
   "De Vleeschouwer 2009 Part II's GLUTAMINE column" was seen, because
   `k1_kinetic_parameters.md` sec. 2c "prints Gln and Cys side by side" and
   `k3_final_parameter_inventory.md` sec. B5.5 prints the glutamine promotion percentages. **The
   exposure is therefore not new in kind, but this dossier raises its fidelity**, and the honest
   handling is to say so: no glutamine number here may enter a parameter, a bound, an
   initialisation or a fit row, the literal-grep firewall in
   `tests/unit/test_kinetic_core_b3.py` covers executable code and not this file, and **this
   dossier must not be cited as a `dossier_anchor` for anything on the glutamine side.** The
   already-recorded observation that glutamine's promotion grows with temperature in Claeys'
   liquid pot (154.9 -> 321.7 %) and shrinks at a_w 0.92 (267.2 -> 120.0 %) is confirmed here from
   the primary tables; **it is not a new finding of this dossier** and it changes nothing about the
   hold-out.
2. **`devleeschouwer2009b.pdf` is a DIFFERENT PAPER and it has no dossier.** It is
   **"...Part I: Effect of the type of sugar", Food Chemistry 114 (2009) 116-126**, same four
   authors, same submission and revision dates, doi 10.1016/j.foodchem.2008.09.024, cited
   throughout Part II as "in press". Its abstract studies glucose, fructose and sucrose in equimolar
   asparagine-sugar systems at the same a_w 0.92 and the same 120-200 C. **This matters because
   three shipped registry rows — `k_asn_glc`, `k_int1_acr`, `k_asn_asp` — come from Part I's
   Table 3 glucose column, and Part I's FRUCTOSE and SUCROSE columns are declared hold-outs.** No
   dossier is written for it here, as instructed. Anyone who writes one must handle the two
   hold-out columns the way Flags 1 handles the glutamine column.
3. **The registry's citation for this paper is wrong in two places.** `_DV2_SOURCE` reads
   "De Vleeschouwer et al. 2009 Part II, **J Agric Food Chem 57:539-546**, Table 3 p. 542". The
   correct citation is **Food Chemistry 114 (2009) 535-546** — wrong journal, wrong volume, wrong
   start page. **Table 3 is on p. 542, so that part is right.** The `dossier_anchor` strings point
   at inventory sections and are unaffected. Fix the four `_DV2_SOURCE` rows.
4. **The three systems are not at the same concentration, and half of Table 3 was transferred
   across that difference.** Table 1's moistures imply 2.74, 1.97 and 1.85 mol/kg for the control,
   glutamine and cysteine systems (mine, section 3 arithmetic 2). The ODE integration uses each
   sample's *measured* initial concentrations, so that is handled; but **k_INT (second order),
   k_B (second order), k_M, k_C and k_Asp were fixed at estimates made at 1.4x the concentration**,
   and the same authors' own De Vleeschouwer 2008b studied "different initial reactant
   concentrations and ratios" precisely because that can matter. The consistency conclusion in
   section 3.3.2 is therefore a test of the *model structure* at a fixed parameter set, not a
   demonstration that the constants are concentration-independent. The registry's
   `_DV2_CONDITIONS` string says "~3 mol/L"; **for the cysteine system 1.85 mol/kg is the better
   figure**, and it is a factor of 1.6 on every bimolecular flux derived from it.
5. **No initial concentration is printed anywhere in this paper**, and the hydration state of the
   asparagine and the glucose is not stated. Every concentration is in Figs. 1 and 2 only. My
   derivation in section 3 arithmetic 2 assumes the Karl-Fischer moisture is w/w of the equilibrated
   powder, a density of 1 kg/L, and anhydrous molecular weights; the anhydrous-vs-monohydrate choice
   moves the absolute numbers by 5-6 % and the ratios by less than 2 %.
6. **The melanoidin response is wrong by the authors' own account, in both columns and in opposite
   directions.** Browning is converted with the asparagine-glucose eps of 282 L/mol cm in a system
   that also makes glutamine-glucose melanoidins (eps 498, i.e. **1.77x higher**) or
   cysteine-glucose melanoidins (no eps exists; Ashoor & Zent put the colour at ~25 % of
   asparagine-glucose), plus caramel. The model **underestimates** browning for glutamine and
   **strongly overestimates** it for cysteine. The paper's own defence is the important part and
   should travel with any use of the table: recalculating with the glutamine eps changed **none** of
   the estimated parameters significantly, "which proves that the kinetic parameters describing
   acrylamide formation and elimination reactions are independent of the extinction coefficient
   used". That defence covers `k_F`, `k_E` and `k_E2`; it does **not** cover `k_INT2`, which is
   fitted directly against the browning response and which the authors themselves mark as apparent
   and print in grey.
7. **Table 3 carries an undefined footnote marker.** Ea_F in the cysteine column is printed as
   **159.2\*** while every other fixed parameter in that column carries footnote **c**. The table
   defines footnotes a, b and c only; there is no `*` anywhere else in the paper except the
   corresponding-author mark. The likeliest readings are (i) a typo for `c`, or (ii) a marker for
   the section 3.3.2 observation that the *freed* re-fit gave a higher Ea_F than the control. **The
   registry ships Ea for `k_int1_acr` as 159.2 ± 29.5 from Part I, so nothing depends on resolving
   this**, but it should be resolved before the number is quoted from Part II.
8. **Two printed typos in the glutamine ODE set (Eqs. 2-12, p. 541), both caught by comparison with
   the cysteine set (Eqs. 13-23, p. 543).** (i) **Eq. (8) reads `d[AA]/dt = k_B·[Int1] - k_E·[AA]`,
   where it must be `k_F·[Int1]`** — Eq. (19) has `k_F` in the same position, `k_B` is the browning
   constant with a different unit (M^-1 min^-1 against min^-1), and Scheme 1 puts `k_F` on the
   Int1 -> AA arrow. (ii) **Eq. (7) reads `d[Browning]/dt = k_B[Int2][Asn] + k_Int2·[Glc]·[Glc]`,
   where it must be `k_INT2·[Gln]·[Glc]`** — Eq. (18) has `k_INT2[Cys][Glc]`, and Eq. (3) already
   uses `k_INT2[Gln][Glc]` in the glucose balance. Neither typo affects the fitted values, which
   come from the software's own model; both would corrupt any re-implementation typed from the
   paper. **Anyone rebuilding this network from the printed equations must use the cysteine set as
   the template.**
9. **Table 2's cysteine column has one significant figure and its glutamine column repeats a value.**
   0.2 / 0.2 / 0.2 / 0.1 / 0.1 cannot resolve a temperature trend, so it cannot be used to test
   Ea_E2 = 51.3 against Ea_E = 113.2. The glutamine column prints **132.4 at both 160 and 180 C**,
   identical to four figures across a 20 K interval in an otherwise monotone series
   (267.2, 180.1, 132.4, 132.4, 120.0) — possibly real, possibly a duplicated cell.
10. **What this paper does NOT contain**: any absolute acrylamide concentration or any other
    concentration (all figure-only); any pH; any second water activity; any competitor other than
    glutamine and cysteine; any identification of the acrylamide degradation products `DP`, of the
    aspartic-acid sink `X`, or of the cysteine sink `Y`; any uncertainty on any fixed parameter; any
    detection limit for acrylamide; any test of the assumed second order of the scavenging step; any
    supplementary material.
11. **What to request from the authors**: (i) the numeric data behind Figs. 1 and 2 — five
    temperatures x six responses x two systems, which would turn this paper from a constants source
    into a benchmark source, and which is the single largest gap; (ii) the initial molar
    concentrations of each system, which would replace my derivation in section 3 arithmetic 2;
    (iii) the meaning of the `*` on Ea_F; (iv) confirmation of the two ODE typos; (v) the identity
    of `Y`, given that the fitted k_Y = 0.35 min^-1 makes it the dominant cysteine sink at 160 C
    apart from the sugar reaction; (vi) a cysteine concentration series at a_w 0.92, which is the
    experiment that would settle the assumed second order inside one laboratory and one matrix.
12. **Registry gaps against `data/keys/compounds.yml`**: `acrylamide` is present and is the only
    compound in this paper that is. **Absent: asparagine, glutamine, cysteine, cystine, aspartic
    acid, glutamic acid, glucose, melanoidins, and the Michael adduct S-(2-carbamoylethyl)cysteine**
    — the last of which is the named product of the registry's `k_acr_cys` transformation string and
    has no id, so the reaction's product cannot be referred to outside the network's own local
    names. The lane's internal species names (`ACR`, `Asn`, `Glc`, `Int1`, `Int2`) are network-local
    and are not registry ids. A benchmark built from this paper would need at least asparagine,
    glucose and cysteine keyed, and cysteine would need an alias policy because **this paper's
    `[Cys]` is cysteine + cystine measured as one number.**
