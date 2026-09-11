# Lee 2022 — EXTRACTION (a real low-moisture sponge-cake matrix carrying ONLY glucose, or glucose + leucine, as declared precursors; baked at oven set-points 140, 170 and 200 C at two fan frequencies for 0-120 min, with thermocouples at three depths, and 12 markers quantified against time: glucose, fructose, free NH2, glucosone, 1-DG, 3-DG, 3,4-DG, glyoxal, methylglyoxal, diacetyl, furfural, HMF and A420 browning)

### THE ONE PAPER IN THIS CLUSTER WITH A CONTROLLED NITROGEN-FREE CONTROL: the glucose-only cake browns with **no amine in the system at all**, and the glucose + leucine cake browns faster and darker — so the pair is a within-study contrast between a brown polymer whose nitrogen content is zero by construction and one that has nitrogen, which is exactly the `MEL_C` / `MEL_N` distinction the trunk draws. It also reproduces Kocadagli & Gokmen's 3-DG : glucosone : 1-DG ratio in a second laboratory and a second matrix.

**Source on disk:** `data/articles/lee2022.pdf` (28 pp. including the HAL cover sheet; the article
is Food Chemistry 2022, **376**, 131917).
Read from the `pdftotext -layout` text layer
(`scratchpad/articles/lee2022.txt`), which is clean throughout — this is the **HAL author
accepted manuscript** (hal-03819621, CC BY-NC 4.0), not the typeset Elsevier version, so its
figures are separated from their captions but its text is complete and its one table came through
intact. **The paper has exactly one table (Table 1, an analytical method table) and no other.**
Figures 1, 3, 4 and 5 carry **every concentration-time datum in the paper** and are
**figure_only**; the levels the authors quote in the running text are carried below and are
sourced to the text, not to the figures. **The underlying kinetic data are deposited under six
separate DOIs on the INRAE data portal** and are not on disk — see Flags 1, which is the most
actionable flag in this dossier. Repo status before this dossier: `lee2022.pdf` has **no
extraction dossier** and is not cited anywhere in `src/`.

**A note on the file names.** `lee2022.pdf` and `lee2024.pdf` share the first half of their title
and were checked against each other before this dossier was written. **They are two different
papers**, same authors, same model cake, same baking rig, different journals and different
measured species: this one measures the non-volatile precursors, alpha-dicarbonyls and furanic
compounds; `lee2024.pdf` (Food Research International 2024, 183, 114183) measures the **volatile**
markers extracted on-line during baking. Both dossiers exist.

## 0. Identity

| field | value |
|---|---|
| Title | "Unravelling caramelization and Maillard reactions in glucose and glucose + leucine model cakes: Formation and degradation kinetics of **precursors, alpha-dicarbonyl intermediates and furanic compounds** during baking" |
| Authors | J. (Jeehyun) Lee; **S. (Stéphanie) Roux — corresponding**; E. (Even) Le Roux; S. (Séverine) Keller; B. (Barbara) Rega; C. (Catherine) Bonazzi. Université Paris-Saclay, INRAE, AgroParisTech, **UMR SayFood**, Massy / Palaiseau, France |
| Venue | Food Chemistry 2022, 376, 131917 |
| DOI | 10.1016/j.foodchem.2021.131917. Preprint: HAL hal-03819621 |
| Naming | **G** = the glucose-only model cake (a caramelisation-like system, **no nitrogen source**). **G+L** = glucose + leucine. Markers are abbreviated **GCO** (glucosone), **1-D** (1-deoxyglucosone), **3-D** (3-deoxyglucosone), **3,4-D** (3,4-dideoxyglucosone), **GO** (glyoxal), **MGO** (methylglyoxal), **DA** (diacetyl), **F** (furfural), **5-HMF**, **3-MB** (3-methylbutanal). Concentrations are **per gram of dry matter** (`gDM`) throughout |
| Lineage | the **SayFood model-cake programme**: Fehaili 2010 (the instrumented pilot oven), Lee 2020 (the model cake itself and the free-NH2 titration), Srivastava 2018 (the same cake, browning localisation), Cepeda-Vázquez 2018, Ait Ameur 2008. Its kinetic frame is taken from **Kocadagli & Gokmen 2016** (the glucose/wheat-flour multiresponse model), Martins & van Boekel 2003/2005, Brands & van Boekel 2001/2003, De Vleeschouwer 2009 and Goncuoglu Tas & Gokmen 2017 — **five sources the repository already holds dossiers for** |
| Companion on disk | **`lee2024.pdf` is the volatile-marker sequel** — dossier `lee2024_extraction.md` |
| Companions on disk | `kocadagli2016foodchem_extraction.md` (Food Chemistry 211:892 — the paper this one reproduces a ratio from, and one of the B7 furanic block's two sources), `kocadagli2016jafc_extraction.md`, `martins2005_extraction.md` and `martins2003*_extraction.md` (the trunk's own network), `brands2001_extraction.md`, `devleeschouwer2009_extraction.md`, `goncuoglutas2017_extraction.md`, `nguyen2016_extraction.md`, `hollnagel2000_extraction.md` |

## 1. Why it matters

**(a) It is the cluster's only nitrogen-free browning control, and that is what it contributes to
the C/N question.** The trunk carries the brown polymer as two elemental pools and
`melanoidin_c_over_n` in `src/kinetic_core/species.py` returns **NaN when `MEL_N` is zero or
below**, precisely because a carbon-only polymer has no C/N. The **G model cake is the physical
realisation of that state**: a solid food matrix with glucose as the *only* declared reactive
precursor, no amino acid, no protein — the starch, methylcellulose and HPMC are the structure —
and it browns measurably. The paper says the consequence in as many words: "**Melanoidins differ
from the brown polymers previously described in that they contain nitrogen**", and
"the difference in browning levels between G and G+L also indicated ... different compositions
and characteristics for brown compounds (brown polymers from caramelization or melanoidins from
Maillard reaction)". So the G / G+L pair is a **within-study contrast at the two extremes of the
diagnostic**: one polymer with C/N undefined, one with C/N finite, both measured for browning on
the same instrument in the same matrix on the same day. Neither is measured elementally
(Flags 5), so this cannot become a number — but it is the correct physical picture of what the
trunk's two-pool bookkeeping is for.

**(b) It carries an amine the trunk cannot substitute, and the reason is arithmetic.**
`species.py`'s note on `Gly` says an alpha-amino acid "with the same carbon/nitrogen count
substitutes without changing the bookkeeping, but NOT without changing the rates". **Leucine is
C6H13NO2 — six carbons to one nitrogen, against glycine's two.** Under the trunk's own step-9
stoichiometry (`MELANOIDIN_REPEAT_UNIT_CARBON = 8` = 6 from 3-deoxyglucosone + 2 from glycine), a
leucine repeat unit would be **6 + 6 = 12 carbons per nitrogen** if the amine arrives intact, or
6 + 5 = 11 if it arrives decarboxylated. **So the C/N of a melanoidin is not a property of
"the Maillard reaction"; it is set by which amine is in the pot.** This is the sharpest constraint
in the cluster on how far Mundt 2004's measured 7.64 can be generalised: it is a *glycine* number.

**12 is only the upper bound, and the volatile companion measures the other end.** If leucine is
fully Strecker-degraded it contributes **no carbon at all** to the polymer — the carboxyl leaves
as CO2, the side chain plus alpha carbon leave as the C5 Strecker aldehyde **3-methylbutanal**,
and only the nitrogen transfers to the dicarbonyl — putting the repeat unit at **6 + 0 = 6 C per
N**. `lee2024.pdf` finds 3-methylbutanal is "**the most abundant of all the reaction markers**",
escaping into the oven vapour and off-scale at 200 C. **The honest bracket for a leucine
melanoidin is therefore 6 to 12, and the measured evidence favours the low end.** What decides
where in that bracket a real amine sits is the **volatility of its Strecker aldehyde**: glycine's
is formaldehyde, which Mundt's radiochemistry shows is *retained* in the polymer, whereas
leucine's is measured leaving. The three-row bracket for both amines is worked in
`lee2024_extraction.md`, section 3 arithmetic 4.

**(c) It reproduces a ratio the repository already uses, in a second laboratory and a second
matrix.** The paper reports 3-D formed **25-fold** more than GCO and **10-fold** more than 1-D,
and says so explicitly against Kocadagli & Gokmen 2016, who "reported 3-D concentrations that
were roughly **26-** and **9-fold** higher than those of GCO and 1-D, respectively, in a low
moisture wheat flour/glucose model system". Kocadagli & Gokmen 2016 is one of the two sources of
the trunk's B7 furanic block (`src/kinetic_core/parameters_furanic.py`, and the `INT` species
note in `species.py` records its constraints). **A cross-laboratory, cross-matrix reproduction of
a dicarbonyl ratio the repository's parameters depend on is a scorecard-relevant finding**, and it
comes with a real-food matrix, which the scorecard's stated gap asks for.

**(d) It is a low-moisture, real-matrix, non-isothermal system — which is what the trunk lacks
and also what makes it unusable as a benchmark.** The corpus is dominated by stirred aqueous pots
and freeze-dried powders. This is a 1.3 cm sponge cake with a measured internal temperature
gradient, a moisture gradient, an alveolar structure and a crust. The paper measures the
temperature at three depths every 2 s. That is exactly the matrix the repository's data wishlist
asks for **and** exactly the reason nothing here can be fitted without solving heat and mass
transport first — which the authors themselves name as future work: "It will require specific
modeling work by coupling transfer equations and chemical reaction models."

**(e) It supplies a caveat on the browning readout that applies directly to the trunk's
hold-out.** `results/validation/kinetic_core_b1_holdout_report.json` scores the melanoidin
trajectory against Martins' A470 with eps = 0.64 L/(mmol·cm). Lee measures A420 on a
**water-soluble** extract and observes that "at 200 C, the measured absorbance **decreased** at
the end of baking, but this was probably due to high molecular weight and **water-insoluble**
melanoidins and brown polymers **that were not covered by the measurement**". A soluble-extract
absorbance is not a measure of total polymer once the polymer starts precipitating — the same
failure Knol 2005 records for his A470 at 160-200 C. **Any browning hold-out has an upper
severity beyond which the readout under-reports, and this paper locates it in a real matrix.**

**(f) It prints two structural conclusions that a network can act on.** From the timing and
magnitudes at 170 C the authors conclude: **the 1-D -> glyoxal and 1-D -> diacetyl pathway "could
be set aside as it was negligible"**; **furfural comes predominantly via 3-D, not via fructose**;
and **5-HMF comes mainly from fructose, with 3,4-D contributing only "a minor addition"**. All
three are topology statements about limbs the trunk's B7 block carries.

What this paper does NOT give the repository: any rate constant, any activation energy, any
tabulated concentration, any elemental analysis, any C/N, any measurement of 3-methylbutanal
(that is the 2024 companion), and any isothermal condition.

## 2. Methods as they matter to a model

- **The matrix.** A model sponge cake built from purified water, **native corn starch (12.4 % w/w
  water**, Cargill), food-grade **methylcellulose (MC, type SGA7C)** and **hydroxypropyl
  methylcellulose (HPMC, type K250M)** (Dow). A hydrocolloid solution of HPMC and MC is foamed
  with corn starch and enriched **immediately before baking** with the declared precursors. The
  cake "imitates sponge cake in terms of its alveolar structure and manufacturing method (using
  same operations of mixing, foaming and baking)". **The point of the design is that glucose, or
  glucose plus leucine, are the only reactive precursors**; the paper never states their
  concentrations in the formula directly, but reports "an initial glucose concentration close to
  **2 mmol·gDM^-1**".
- **The two formulations.** **G** = glucose only. **G+L** = glucose + L-leucine. Glucose from
  Roquette Frères, leucine from Sigma-Aldrich, both food grade, >= 99 %. **The leucine
  concentration is never printed anywhere in this paper** (Flags 4).
- **Geometry and load.** 20 g of batter per disposable cylindrical aluminium mould, diameter
  **6.6 cm**, mould height 4 cm, **cake height 1.3 cm**. **14 cakes** baked together at 140 C,
  **16** at 170 or 200 C.
- **Baking.** An **instrumented pilot oven with precise and uniform temperature control**
  (Bongard, Wolfisheim; the rig of Fehaili 2010). Three oven set-points — **140, 170, 200 C** —
  crossed with two **fan frequencies, 25 and 50 Hz**. **The 170 C trials were run in
  triplicate; 140 C and 200 C were n = 1 per fan setting** (Flags 3).
- **Sampling.** 0, 6, 12, 24, 37, 56, 75 and 90 min. At 140 C two extra times (105, 120 min) for
  G; for G+L at 140 C the 90 min point was replaced by 95 min and 120 min added. The long
  durations were chosen "in order to push the degradation of the substrates and to reveal the
  behaviors of the intermediates as completely as possible".
- **Quenching.** Sealed in tared aluminium bags immediately, plunged into a water/ethanol bath at
  **−8 C for at least 15 min**, stored at −20 C for 24-48 h, frozen at −80 C for 24 h,
  freeze-dried (ice condenser −85 C), ground 1 min, split into eight aliquots, stored at −20 C.
- **The temperature history, which is the crux.** Three **1 mm J-type thermocouples**,
  individually calibrated against a certified probe, logged every **2 s** by LabVIEW, at three
  positions in one cake: **bottom** (central axis, cake/mould interface), **centre** (central
  axis, half height) and **surface** (periphery, at the cake/air interface). Findings printed in
  the text: the surface rises asymptotically towards the oven set-point; the centre rises slowest
  and shows **a plateau whose value and duration depend on the set-point**, corresponding to the
  balance between energy supplied and energy consumed by water evaporation; the bottom lies
  between. **The temperature becomes uniform and equal to the oven only after 120 min at 140 C,
  90 min at 170 C and 70 min at 200 C.** Times to reach 120 C — the threshold the paper quotes
  for glucose caramelisation — are **surface / centre: 5 / 35 min at 200 C, 8 / 55 min at 170 C,
  20 / 80 min at 140 C**.
- **A negative result on the two design variables.** "**No significant differences between the
  two formulations and convection levels (fan frequencies) could be detected on the temperature
  profiles nor on the water content profiles**", so the leucine addition changed no heat or mass
  transfer, and the fan frequency changed neither. In the marker figures the two fan settings are
  therefore treated as repetitions of one temperature.
- **Moisture.** Dry matter by desiccation 24 h at 105 C; cake moisture from the mass difference
  between the initial dough (10^-2 g) and the sampled cake (10^-4 g). **No moisture or
  water-activity value is printed anywhere** (Flags 6).
- **Browning.** Water-soluble brown pigments extracted at 20 C, **0.8 g freeze-dried sample per
  25 g ultrapure water**, vortex 30 s at 3000 rpm, filtered at 8 um (Whatman), **absorbance at
  420 nm** (Secomam S250), diluted to keep A < 1, then corrected for dry matter and dilution.
  **Water-soluble only** — see (e) above and Flags 7.
- **Glucose and fructose.** Extracted at 20 C (0.4 g / 20 g water), vortex, 5000g 15 min, Carrez I
  and II clarification, 10000g 2 min, 0.2 um nylon. **UHPLC-CAD** (Dionex U3000 + Corona Veo RS
  Charged Aerosol Detector), Acquity BEH Amide 100 x 2.1 mm 1.7 um at 35 C, water/acetonitrile
  both with **10 mM NH4OH**, 0.26 mL/min, 19 min run, 2 uL injected. External standard, 4-6
  points from 0.01 to 3 g/L, **quadratic** response, no forced zero.
  **LOD/LOQ: glucose 7.803 / 26.01 mg/L; fructose 9.339 / 31.13 mg/L.**
- **Free NH2 (the leucine readout).** **Titration against 0.005 N NaOH** on a Metrohm 809
  Titrando, per Lee 2020; expressed as **moles of NH2 per gram of dry matter**. Linearity checked
  at six points over 0.0039-0.0382 mmol/g. **LOD/LOQ 0.761 / 2.536 umol/g.** Note that this is a
  free-amino-group titration, **not** a leucine assay (Flags 4).
- **alpha-Dicarbonyls.** Extracted at 20 C (0.4 g / 20 g water), 5000g 15 min, then **5 mL
  supernatant + 5 mL acetonitrile** and 7000g 3 min "**in order to precipitate colloidal
  polymers**". Derivatised: 900 uL supernatant + 150 uL phosphate buffer pH 7.0 (100 mmol/L) +
  150 uL **o-phenylenediamine 0.2 % in DETAPAC (10 mmol/L)**, **60 C for 30 min**, cooled 1 min,
  then **100 uL of quinoxaline-5,6,7,8-d4 (10 mg/L) as internal standard**, filtered 0.2 um.
  **UHPLC-MS-QToF** (Waters Acquity H-Class + Xevo G2-XS QTof, ESI+), Purospher STAR RP-18e
  150 x 4.6 mm 5 um at 30 C, water/methanol both with 1 % formic acid, 0.3 mL/min, 15 min run,
  1 uL injected; resolution mode, 0.5 s scan, 50-600 m/z centroid, internal lock-mass on leucine
  enkephalin. **Note the internal standard is deuterated and added AFTER derivatisation**, so it
  corrects for injection and ionisation but not for derivatisation yield.
- **The quantification shortcut that governs four of the seven dicarbonyls.** Calibration at 8
  points against authentic **quinoxaline, methylquinoxaline and dimethylquinoxaline** only.
  "**The calibration curve for quinoxaline was used to quantify several commercially
  non-available compounds (glucosone, 1-deoxyglucosone, 3-deoxyglucosone and 3,4-dideoxyglucosone
  derivatives).**" So **GCO, 1-D, 3-D and 3,4-D are all quantified against the response factor of
  unsubstituted quinoxaline**, and their absolute scale carries an unknown multiplicative error —
  the same author-declared semi-quantitation the trunk already records for `DDG` from the
  Kocadagli papers (`species.py`, the `DDG` note, K5a C22). See Flags 2.
- **Furfural and 5-HMF.** From the same extract as the sugars, stored at most 2 weeks at −20 C.
  Thermo U3000 + PDA, Acquity HSS T3 100 x 2.1 mm 1.8 um at 40 C, water/acetonitrile both with
  0.1 % formic acid, 0.5 mL/min, 11 min run, 1.4 uL injected; **detection at 277 nm (furfural) and
  284 nm (5-HMF)**. External standard, 7 points 0.05-21.1 mg/L, linear, no forced zero.
  **LOD/LOQ: furfural 5.93 / 19.8 ug/L; 5-HMF 0.190 / 0.634 ug/L.**
- **Repeatability.** For 3-D, **eight successive extractions from the same freeze-dried sample,
  RSD = 0.01 %** — an extraction-repeatability figure, not a between-bake figure (Flags 3).

## 3. Tables re-typed

### Table 1. "Derivative forms of alpha-dicarbonyl compounds associated with their retention times and selected quantifier ions"

This is **the only table in the paper**.

| compound | retention time (min) | selected ion (m/z) |
|---|---:|---:|
| Glucosone (GCO) derivative | 3.52 | 235.1077 |
| 1-Deoxyglucosone (1-D) derivative | 4.10 | 217.0971 |
| 3-Deoxyglucosone (3-D) derivative | 4.48 | 217.0971 |
| 3,4-Dideoxyglucosone (3,4-D) derivative | 6.56 | 251.1009 |
| Quinoxaline (glyoxal (GO) derivative) | 7.65 | 131.0598 |
| Methylquinoxaline (methylglyoxal (MGO) derivative) | 9.40 | 145.0760 |
| Dimethylquinoxaline (diacetyl (DA) derivative) | 10.82 | 159.0922 |
| **Internal standard** | | |
| Quinoxaline-5,6,7,8-d4 | 7.54 | 135.0855 |

Note that **1-D and 3-D share the quantifier ion 217.0971** and are separated only by retention
time (4.10 vs 4.48 min), a 0.38 min gap on a 15 min gradient.

### The PubChem identifiers the paper prints for the compounds it studies

| compound | PubChem CID | IUPAC name as printed |
|---|---|---|
| Glucose | 107526 | — |
| Fructose | 2723872 | — |
| Leucine | 6106 | — |
| Glyoxal | 7860 | oxaldehyde |
| Methylglyoxal | 880 | 2-oxopropanal |
| Diacetyl | 650 | butane-2,3-dione |
| Glucosone | 159630 | (4S,5R)-4,5,6-trihydroxy-2-oxohexanal |
| 1-Deoxyglucosone | 11228966 | (4R,5R)-4,5,6-trihydroxyhexane-2,3-dione |
| 3-Deoxyglucosone | **114839** | (4S,5R)-4,5,6-trihydroxy-2-oxohexanal |
| 3,4-Dideoxyglucosone | 132520491 | (5R)-5,6-dihydroxy-2-oxohexanal |
| Furfural (printed "Furfual") | 7362 | furan-2-carbaldehyde |
| 5-Hydroxymethylfurfural | 237332 | 5-(hydroxymethyl)furan-2-carbaldehyde |

**The IUPAC names printed for glucosone and 3-deoxyglucosone are identical** — "(4S,5R)-4,5,6-
trihydroxy-2-oxohexanal" — which cannot be right for two different compounds. Recorded as
printed; see Flags 9.

### The deposited datasets (DOIs printed in the text, files NOT on disk)

| dataset | DOI |
|---|---|
| Figure 1 — temperature profiles | doi.org/10.15454/MI4A3S |
| water content profiles ("data not shown") | doi.org/10.15454/18HVVK |
| Figure 3 — precursor consumption | doi.org/10.15454/GUXKYQ |
| fructose kinetics | doi.org/10.15454/95JTQS |
| furanic compounds | doi.org/10.15454/IICJV3 |
| browning | doi.org/10.15454/QYCBIW |
| Figure 4 — alpha-dicarbonyl intermediates | doi.org/10.15454/HV3UTK |

### Every concentration and ratio printed in the running text

**All concentration-time data are in Figures 1, 3, 4 and 5 and are figure_only.** The following
are the levels and ratios the authors state in words; they are printed text and are carried as
such.

| quantity | value | conditions | where |
|---|---|---|---|
| initial glucose | **close to 2 mmol·gDM^-1** | both formulations | Results, "fructose" |
| **maximum fructose accumulation** | **did not exceed 0.08 mmol·gDM^-1** | any condition | Results, "fructose" |
| 5-HMF | **20 (G) and 50 (G+L) umol·gDM^-1** | after 80 min at 200 C | Results, "Furanic and brown products" |
| furfural | **2 (G) and 8 (G+L) umol·gDM^-1** | after 80 min at 200 C | same |
| glyoxal (GO) | **not detected at all in G**; **0.05-0.2 umol·gDM^-1** in G+L, "only under the most extreme conditions (200 C, high convection level)", data not shown | | Results, "alpha-dicarbonyl intermediates" |
| glucosone (GCO) in G at 140 C | **not measured** (below detection) | | same |
| glucosone in G at 170 C | rose gradually to **0.03 umol·gDM^-1** after 90 min | | same |
| glucosone in G at 200 C | bell-shaped, **peak at 36 min** | | same |
| **3-D vs GCO** | **3-D 25-fold larger** | this study | Results |
| **3-D vs 1-D** | **3-D 10-fold larger** | this study | Results |
| the same ratios in Kocadagli & Gokmen 2016 | **26-fold** and **9-fold** | low-moisture wheat flour / glucose | Results (**not measured here**) |
| **3-D, G+L vs G** | **roughly 11-fold higher in G+L** | | Results |
| **GCO, 1-D and 3-D, G+L vs G** | **about 4- to 10-fold higher in G+L** | | Results, conclusion of the dicarbonyl section |
| 3,4-D vs 3-D, in G | **same order of magnitude** | | Results |
| 3,4-D, G+L vs G | **same order of magnitude**, but all G+L kinetics bell-shaped with an earlier, lower maximum as temperature rose | | Results |
| MGO, G+L vs G | "very similar", "slightly higher in the presence of leucine" | | Results |
| diacetyl (DA) in G | "only very small amounts ... at 200 C under high convection" | | Results |
| all alpha-dicarbonyls | **order of umol·gDM^-1** against precursors at **order of mmol·gDM^-1** | | Results |
| glucose consumption, G at 200 C | "almost totally consumed after 90 min" | | Results |
| glucose consumption, G+L at 200 C | "complete consumption in less than 80 min" | | Results |
| glucose consumption, G+L vs G at 140 and 170 C | "higher and more rapid" in G+L | | Results |
| leucine (free NH2) consumption | confirmed at all three set-points; kinetics "very similar at 140 C and 170 C", slightly faster at 200 C; **"leucine degradation was significantly lower than that of glucose"** | | Results |
| fructose, G vs G+L at 170 C | **G: bell-shaped** (accumulates then consumed); **G+L: still rising at the end of baking**; and less fructose in G+L up to 56 min | | Results |
| fructose vs temperature | higher concentrations at **170 C than at 200 C** for both formulations | | Results |
| lag phase for the furanics | **~20 min**, similar to fructose | | Results |
| 3,4-D lag behind 3-D | **at least 12 min** | 170 C | Results |
| browning | more pronounced and more rapid in G+L **at all temperatures**; at 200 C the measured A420 **decreased at the end of baking** | | Results |
| times to reach 120 C, surface / centre | **5 / 35 min (200 C), 8 / 55 min (170 C), 20 / 80 min (140 C)** | | Results, from Figure 1 |
| time to a uniform product temperature | **120 min (140 C), 90 min (170 C), 70 min (200 C)** | | Results |

### The topology conclusions the paper draws from the timings (printed, and load-bearing)

- **Glucose -> 1,2-enediol dominates over glucose -> glucosone.** "The pathway from glucose to
  1,2-enediol and subsequently into fructose and 3-D predominated over the pathway to GCO", and
  the rate constant of glucose -> 1,2-enediol "was certainly the highest".
- **Fructose forms faster than 3-D from the shared 1,2-enediol.** Both accumulate from 24 min,
  "but with a noticeably higher quantity of fructose".
- **5-HMF is mainly from fructose; 3,4-D contributes only a minor addition.** "The contribution
  from 3,4-D to 5-HMF probably existed but only made a minor addition to what was formed from
  fructose."
- **Furfural is mainly from 3-D, not from fructose.** "In view of the time lag for its
  accumulation, this pathway [from fructose] could be set aside in favor of that via 3-D."
- **The 1-D -> glyoxal and 1-D -> diacetyl limb is negligible.** "No GO was detected and hardly
  any DA. This means that the degradation pathway from 1-D to GO and DA could be set aside."
- **1-D -> MGO is very fast and MGO -> brown pigments is slow.** "MGO appeared simultaneously
  with 1-D and was even detected in larger quantities than 1-D towards the end of baking. It can
  thus be assumed that the conversion of 1-D to MGO was certainly very rapid and that the
  degradation of MGO towards brown pigments was probably slow."
- **No 2,3-enediol accumulation.** "1-D displayed the same time lag as 5-HMF ... indicating no
  accumulation of 2,3-enediol, as otherwise 1-D would have appeared later."
- **Leucine adds pathways rather than changing constants.** "The rate constants of the reactions
  described above remain unchanged when leucine is added to the system. Only new reaction pathways
  that include leucine are supposed to be added." Everything the paper concludes about G+L follows
  from that assumption (Flags 8).

### Arithmetic on the printed numbers (all mine)

**1. The dicarbonyl pool is a rounding error on the carbon balance.** Initial glucose is
~2 mmol·gDM^-1 = 2000 umol·gDM^-1, i.e. **12 000 umol of carbon per gram DM**. The largest
furanic, 5-HMF at 50 umol·gDM^-1 in G+L, is **300 umol C**, or **2.5 %** of the initial glucose
carbon; every alpha-dicarbonyl is "in the order of umol·gDM^-1", so the whole dicarbonyl pool is
well under **1 %**. **Over 95 % of the glucose carbon in a fully baked cake is unaccounted for by
any measured marker**, and the residue is brown polymer, CO2, water and volatiles. This is the
clearest statement in the corpus of why a pool like `FRAG_C` has to exist.

**2. Fructose never holds much of the flux.** Maximum fructose 0.08 mmol·gDM^-1 against 2
mmol·gDM^-1 initial glucose is **4 %**. Fructose is a fast-turning-over intermediate here, not a
reservoir — consistent with the trunk carrying it as an `intermediate`.

**3. The furanic ratios, G+L versus G (mine).** 5-HMF 50/20 = **2.5x**; furfural 8/2 = **4.0x**.
Both at 80 min, 200 C. So **leucine roughly doubles to quadruples the furanic yield** in this
matrix, and the effect is larger on furfural than on HMF — which fits the paper's own routing
(furfural via 3-D, and 3-D is 11x higher in G+L).

**4. HMF : furfural (mine).** 20/2 = **10 : 1 in G**; 50/8 = **6.25 : 1 in G+L**. The ratio
narrows with leucine, again consistent with furfural's 3-D route being the one leucine boosts.

**5. Does the paper's own 25 : 10 : 1 dicarbonyl ordering hold with the 11x leucine effect?
(mine.)** If 3-D is 25x GCO and 10x 1-D in one formulation, and GCO, 1-D and 3-D are all 4-10x
higher in G+L, then the **ratios among the three are approximately preserved** across
formulations while their absolute levels move together. That is a non-trivial internal
consistency and it is what makes the Kocadagli comparison meaningful: **the 26 : 9 : 1 ordering
survives a change of matrix (wheat flour vs starch/cellulose cake) and a change of amine
(none/leucine vs wheat-flour protein).**

**6. What a leucine melanoidin repeat unit would be under the trunk's own rules (mine).**
`MELANOIDIN_REPEAT_UNIT_CARBON = 8` = 6 (3-DG) + 2 (glycine). Substituting leucine (C6N1) gives
**12 C per N intact, 11 if only decarboxylated, and 6 if fully Strecker-degraded** (in which case
the amine's carbon all leaves — the carboxyl as CO2 and the rest as the C5 aldehyde
3-methylbutanal, which `lee2024.pdf` measures escaping into the oven vapour as the most abundant
of its ten markers). **So the bracket is 6 to 12**, against Mundt & Wedzicha's *measured* glycine
melanoidin C/N of 7.64 — whose corresponding bracket is 8 / 7 / 6, and which lands at 7.64 because
glycine's Strecker aldehyde is formaldehyde and Mundt shows it is retained. Nothing in this paper
measures any of it — this is arithmetic on the trunk's own constant, recorded to show that the C/N
diagnostic is amine-specific, that a single number cannot serve both amines, and that the
deciding variable is the volatility of the amine's Strecker aldehyde. Worked in full in
`lee2024_extraction.md`, section 3 arithmetic 4.

## 4. Numbers the repository can use

**Registry mapping (`data/keys/compounds.yml`, 75 ids).** Of the twelve markers, **three are
keyed**: `hmf`, `furfural` and `2_3_butanedione` (the paper's diacetyl). `3_methylbutanal` is
keyed and is discussed in this paper's reaction scheme as leucine's Strecker aldehyde but is
**not measured here** (it is measured in `lee2024.pdf`). **Glucose, fructose, leucine, glucosone,
1-deoxyglucosone, 3-deoxyglucosone, 3,4-dideoxyglucosone, glyoxal and methylglyoxal are all
absent** from the registry, and four of those (`ODG`, `TDG`, `MGO`, `DDG`) are trunk state
variables.

**Governing conditions on every row: a solid model sponge cake (corn starch + MC + HPMC), 20 g
batter, 1.3 cm thick, baked in an instrumented pilot oven at an OVEN SET-POINT of 140, 170 or
200 C with a fan at 25 or 50 Hz; the PRODUCT temperature is a measured non-isothermal profile with
a surface/centre gradient (see section 2); glucose ~2 mmol·gDM^-1 with or without leucine at an
unstated concentration; 0-120 min; freeze-dried before analysis.**

| quantity | value | unit | conditions | source location | evidence class |
|---|---|---|---|---|---|
| initial glucose | ~2 | mmol·gDM^-1 | both formulations, t = 0 | Results p. 14 | level_only |
| maximum fructose | ≤ 0.08 | mmol·gDM^-1 | any condition | Results p. 14 | level_only |
| 5-HMF | 20 (G) / 50 (G+L) | umol·gDM^-1 | 80 min, 200 C set-point | Results p. 14-15 | level_only |
| furfural | 2 (G) / 8 (G+L) | umol·gDM^-1 | 80 min, 200 C set-point | Results p. 15 | level_only |
| glyoxal | **not detected in G**; 0.05-0.2 in G+L | umol·gDM^-1 | 200 C, 50 Hz only | Results p. 15 | level_only (**a measured null in G**) |
| glucosone | 0.03 | umol·gDM^-1 | G, 170 C, 90 min | Results p. 15 | level_only |
| **3-D : GCO** | **25 : 1** | — | this study, both formulations | Results p. 16 | **within_study_ratio** |
| **3-D : 1-D** | **10 : 1** | — | this study | Results p. 16 | **within_study_ratio** |
| the same in Kocadagli & Gokmen 2016 | 26 : 1 and 9 : 1 | — | wheat flour/glucose, low moisture | Results p. 16 (**cited, not measured**) | level_only (borrowed) |
| **3-D, G+L : G** | **~11 : 1** | — | as above | Results p. 16 | **within_study_ratio** |
| **GCO, 1-D and 3-D, G+L : G** | **4-10 : 1** | — | as above | Results p. 16 | **within_study_ratio** |
| 3,4-D, G+L : G | ~1 : 1 (same order of magnitude) | — | as above | Results p. 16 | within_study_ratio |
| MGO, G+L : G | slightly > 1 | — | as above | Results p. 16 | within_study_ratio (**no number**) |
| 5-HMF, G+L : G | 2.5 : 1 | — | 80 min, 200 C | derived (mine) | within_study_ratio |
| furfural, G+L : G | 4.0 : 1 | — | 80 min, 200 C | derived (mine) | within_study_ratio |
| 5-HMF : furfural | 10 : 1 (G); 6.25 : 1 (G+L) | — | 80 min, 200 C | derived (mine) | within_study_ratio |
| 3,4-D lag behind 3-D | ≥ 12 | min | 170 C | Results p. 18 | within_study_ratio (temporal) |
| furanic lag phase | ~20 | min | all conditions | Results p. 15 | within_study_ratio (temporal) |
| dicarbonyl pool vs precursor pool | umol·gDM^-1 vs mmol·gDM^-1, i.e. **< 1 %** | — | all conditions | Results p. 16; ratio mine | within_study_ratio |
| measured markers as a share of initial glucose carbon | **< 5 %** even at the most severe condition | — | 80 min, 200 C, G+L | derived (mine) | derived_assumption |
| **surface / centre time to 120 C** | 5 / 35 (200 C); 8 / 55 (170 C); 20 / 80 (140 C) | min | measured by thermocouple | Results p. 12, from Figure 1 | **measured level** (a temperature history, not a rate) |
| time to a uniform product temperature | 120 (140 C); 90 (170 C); 70 (200 C) | min | same | Results p. 12 | measured level |
| leucine and fan frequency have **no** effect on the temperature or water profiles | — | — | all conditions | Results p. 12 | level_only (**a measured null; it is what licenses pooling the two fan settings**) |
| **browning is faster and stronger in G+L than in G** | — | A420 on a water-soluble extract | all three set-points | Results p. 15 | level_only (**figure_only for the values**) |
| **A420 falls at the end of baking at 200 C** because the polymer becomes water-insoluble | — | — | 200 C | Results p. 15 | level_only (**a stated readout failure**) |
| topology: 1-D -> GO and 1-D -> DA | **negligible, can be set aside** | — | 170 C, both formulations | Results p. 18; Conclusions | level_only (structural) |
| topology: furfural mainly via 3-D, not fructose | — | — | as above | Results p. 18; Conclusions | level_only (structural) |
| topology: 5-HMF mainly from fructose; 3,4-D a minor addition | — | — | as above | Results p. 18; Conclusions | level_only (structural) |
| topology: no 2,3-enediol accumulation | — | — | as above | Results p. 18 | level_only (structural) |
| topology: 1-D -> MGO fast; MGO -> brown polymer slow | — | — | as above | Results p. 18; Conclusions | level_only (structural) |
| extraction repeatability, 3-D | RSD 0.01 | % over 8 extractions of one sample | — | Methods p. 7 | measured level (**not a between-bake precision**) |
| LOD / LOQ, glucose | 7.803 / 26.01 | mg·L^-1 | UHPLC-CAD | Methods p. 5 | measured level |
| LOD / LOQ, fructose | 9.339 / 31.13 | mg·L^-1 | UHPLC-CAD | Methods p. 5 | measured level |
| LOD / LOQ, furfural | 5.93 / 19.8 | ug·L^-1 | UHPLC-PDA 277 nm | Methods p. 7 | measured level |
| LOD / LOQ, 5-HMF | 0.190 / 0.634 | ug·L^-1 | UHPLC-PDA 284 nm | Methods p. 7 | measured level |
| LOD / LOQ, free NH2 | 0.761 / 2.536 | umol·g^-1 | NaOH titration | Methods p. 6 | measured level |
| retention times and quantifier ions, 7 quinoxalines + IS | see section 3 | min, m/z | UHPLC-MS-QToF ESI+ | Table 1 p. 7 | level_only (method transfer) |
| every concentration-time course, temperature profile and browning curve | — | — | — | Figures 1, 3, 4, 5 | **figure_only** (but deposited — Flags 1) |
| leucine repeat-unit C/N under the trunk's own step-9 rule | **6 to 12** — 12 (intact), 11 (decarboxylated only), 6 (fully Strecker-degraded, amine carbon all lost); `lee2024.pdf` favours the low end | mol C per mol N | — | arithmetic on `species.py` plus Lee 2024's 3-methylbutanal result (mine) | derived_assumption |

### How this bears on the C/N diagnostic and the trunk

**(a) The G cake is a browning system with `MEL_N` identically zero, and it browns.** That is the
one thing this cluster otherwise lacks: every other paper here (Mundt, Fang x2) makes its polymer
from a sugar **and** an amine. Lee's G cake browns from glucose alone in a real low-moisture
matrix. In the trunk's bookkeeping that state has `melanoidin_c_over_n` returning **NaN**, which
is correct and is exactly what `species.py` intends. **The G / G+L pair is therefore the
qualitative validation of the two-pool design**: the same matrix, the same oven, the same
absorbance instrument, and a brown polymer that either has nitrogen or does not.

**(b) But nothing here measures the polymer's composition, so it cannot check the level.** No
CHN, no C/N, no isolation of the pigment, no molecular weight — only A420 on a water-soluble
extract, which the authors themselves say stops tracking the polymer at high severity. **The
comparison object for the trunk's diagnostic remains Mundt & Wedzicha's 7.64 ± 0.21**, not
anything in this paper.

**(c) The amine identity is the finding that most constrains how Mundt's number is used.**
Leucine brings six carbons per nitrogen, glycine two. Any melanoidin C/N is amine-specific, and
the trunk's floor of 8.0 is a *glycine* floor. A future wave that wants a C/N target for a
leucine or a mixed-amine system cannot reuse 7.64.

**(d) Two things could become benchmark rows and one could not.** The **G+L : G ratios** (3-D
11x, GCO/1-D/3-D 4-10x, 5-HMF 2.5x, furfural 4x) are within-study, matrix-matched and
temperature-matched, and they survive the semi-quantitation problem because a response-factor
error cancels between two formulations measured the same way. The **3-D : GCO : 1-D ordering**
(25 : 10 : 1) is a second-laboratory reproduction of a ratio the B7 block's source already
carries. What cannot become a row is **any absolute concentration**, because they are all
figure-only and all rest on the quinoxaline response factor (Flags 2).

**(e) Nothing here is a rate, and the barrier is transport, not data.** The paper's own
conclusion is that turning this into constants "will require specific modeling work by coupling
transfer equations and chemical reaction models to identify all the reaction constants, taking
into account the temperature and water gradients in the products". A cake with a 30-50 min lag
between surface and centre reaching 120 C cannot be fitted by an isothermal integrator. **If the
repository ever wants this dataset, it needs the deposited files (Flags 1) and a transport
model.**

## 5. Flags

1. **The complete kinetic dataset is deposited under seven INRAE DOIs and is NOT on disk.** Every
   figure in this paper carries a `doi.org/10.15454/...` pointer (listed in section 3), including
   one for the water-content profiles that are described as "data not shown". **This is the single
   most valuable follow-up in this dossier**: it would convert every figure_only row above into a
   numeric time series with n = 3 at 170 C, in a real low-moisture matrix, with a measured
   temperature history — which is precisely what `results/validation/data_wishlist.md` asks for.
   Fetch all seven.
2. **Four of the seven alpha-dicarbonyls are quantified against the wrong standard, by the
   authors' own statement.** "The calibration curve for **quinoxaline** was used to quantify
   several commercially non-available compounds (glucosone, 1-deoxyglucosone, 3-deoxyglucosone and
   3,4-dideoxyglucosone derivatives)." Their derivatives are larger, more polar and differently
   ionised than unsubstituted quinoxaline in ESI+. **Every absolute GCO, 1-D, 3-D and 3,4-D level
   in this paper carries an unknown multiplicative error**, and comparisons of those four against
   MGO or DA (which do have authentic standards) are not on one basis. Ratios among the four, and
   ratios of one of them between G and G+L, are unaffected. This is the same defect
   `species.py` already declares for `DDG` from the Kocadagli papers.
3. **The replication is uneven and the quoted precision is the wrong precision.** Only the 170 C
   trials are triplicate; 140 C and 200 C are **n = 1 per fan setting**, and the figures treat the
   two fan settings as repetitions on the strength of the temperature-profile null. The one
   precision figure given — **RSD 0.01 % for 3-D** — is eight extractions of a *single freeze-dried
   sample*, i.e. it measures the extraction and the instrument, **not** bake-to-bake variability.
   No between-bake error is reported anywhere.
4. **The leucine concentration is never printed.** Glucose is given as "close to 2 mmol·gDM^-1";
   leucine has no number anywhere in the paper. Worse, the leucine readout is a **free-amino-group
   titration** expressed as mol NH2 per gram DM, which measures the amine function and not
   leucine, and which the paper itself says reads "the resultant of the amount of consumed leucine
   and that **regenerated** during the Maillard reaction". So the leucine axis has no absolute
   scale and its measured trajectory is a net quantity. **Request the formula from the authors or
   from Lee 2020.**
5. **No composition of the brown polymer is measured — only A420 on a soluble extract.** This
   paper is in the C/N cluster because of the G / G+L contrast and the leucine arithmetic, not
   because it measures anything elemental. It contains no CHN, no isolation, no C/N and no
   molecular weight.
6. **No moisture value, no water activity.** The dry-matter method is described and the water
   profiles are deposited, but **not one moisture or a_w number appears in the text**, and the
   water-content profiles are "data not shown". For a paper whose whole point is a low-moisture
   solid matrix, that is a real gap, and it is the axis on which this matrix differs most from
   the aqueous corpus.
7. **The browning readout fails at the top of the severity range and the authors say so.** A420 on
   a water-soluble extract "decreased at the end of baking" at 200 C because the polymer became
   high-molecular-weight and water-insoluble. **Do not read the late-time browning decline as
   melanoidin destruction.** The same artefact is recorded independently in Knol 2005 (light
   scattering by insoluble particles at 160-200 C).
8. **The whole G+L analysis rests on an unverified assumption stated in one sentence.** "The rate
   constants of the reactions described above remain unchanged when leucine is added to the
   system. Only new reaction pathways that include leucine are supposed to be added." Every
   inference the paper draws about which leucine pathway is active follows from that. It is
   plausible and it is untested here — leucine changes the pH, the ionic environment and the
   water-binding of the matrix, none of which is measured.
9. **The printed IUPAC names for glucosone and 3-deoxyglucosone are identical** —
   "(4S,5R)-4,5,6-trihydroxy-2-oxohexanal" — in the "chemical compounds studied" list. One of the
   two is wrong. Their PubChem CIDs (159630 and 114839) differ, so the identification is not in
   doubt; the printed name is. Also note "Furfual" for furfural in the same list.
10. **1-D and 3-D share a quantifier ion.** Both are extracted at m/z 217.0971 and are separated
    only by 0.38 min of retention time on a 15 min gradient. Any co-elution or peak-shape change
    across the baking series would cross-contaminate the two, and the **3-D : 1-D = 10 : 1 ratio
    is one of this paper's headline numbers**. No chromatogram is shown and no resolution figure
    is given.
11. **"Set-point temperature" is not the reaction temperature and must never be treated as one.**
    140 / 170 / 200 C are **oven air** set-points. The product runs 30-50 min behind at the
    centre, with a measured evaporative plateau, and reaches the oven temperature only at
    120 / 90 / 70 min. Every "at 200 C" in this dossier means "at an oven set-point of 200 C".
12. **What this paper does not contain**: any rate constant; any activation energy; any tabulated
    concentration; any elemental analysis, C/N or polymer composition; any moisture or
    water-activity value; any leucine concentration; any measurement of 3-methylbutanal or any
    other volatile (that is the 2024 companion); any pH; any amine other than leucine; any sugar
    other than glucose; any isothermal condition; any between-bake error estimate.
13. **What to request from the authors**: (i) **the seven deposited datasets** (section 3);
    (ii) the cake formula, especially the leucine and glucose masses per 20 g of batter;
    (iii) the water-content profiles as numbers, and an a_w if one was measured; (iv) the response
    factors, or authentic standards, for the four deoxyglucosone quinoxalines; (v) whether any
    elemental or compositional analysis of the brown pigment was ever done on these cakes.
14. **Registry gaps against `data/keys/compounds.yml`**: `hmf`, `furfural` and `2_3_butanedione`
    are keyed and cover three of the twelve markers. **Glucose, fructose, leucine, glucosone,
    1-deoxyglucosone, 3-deoxyglucosone, 3,4-dideoxyglucosone, glyoxal and methylglyoxal are
    absent** — and four of those are trunk state variables (`ODG`, `TDG`, `DDG`, `MGO`), which is
    the same structural gap every dossier in this cluster records. `3_methylbutanal` is keyed but
    is not measured in this paper.
