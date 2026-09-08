# Pre-registration: wave B18, a pyrazine step on the trunk (written 2026-09-08, not yet run)

## 1. Why

The people this tool is for ask for the roasted note, and no lane makes a pyrazine. The hypothesis
layer reaches the three panel pyrazines from a Strecker charge (rule R28), which says the route
exists; until tonight no dossier on disk held a measured rate for it, so the roadmap held the step
at "rule first, wave when a rate source is read". Five dossiers were read on 2026-09-08
(`zhou2024`, `leahy1989`, `leahy1989a`, `yu2018`, `balagiannis2015`) and the first two carry what a
wave needs: measured formation rates at three temperatures with printed barriers, and within-study
pH ratios from one laboratory.

## 2. The step

On the trunk (glucose + glycine, the only lane whose pot holds glyoxal, methylglyoxal and an
amino acid together), three NET reactions, each the Strecker deamination of the amino acid on the
dicarbonyl followed by the condensation of two aminoketones (rules R07 and R28 of the hypothesis
layer lumped, as Zhou 2024 measured them):

    r_go_pz    2 GO  + 2 Gly -> PZ  + 4 FRAG_C       k_go_pz     rate = k [GO]  [Gly]
    r_mgo_dmp  2 MGO + 2 Gly -> DMP + 4 FRAG_C       k_mgo_dmp   rate = k [MGO] [Gly]
    r_mixed_mp GO + MGO + 2 Gly -> MP + 4 FRAG_C     k_mixed_mp  rate = k [GO]  [MGO]^0 [Gly]  (see below)

PZ is pyrazine, DMP 2,5-dimethylpyrazine, MP methylpyrazine; each glycine leaves its two carbons
as carbon dioxide and formaldehyde, booked to the unassigned fragment pool. The rate law is
Yu 2018's (rate proportional to dicarbonyl times amino acid), second order in L/(mmol·min) at the
trunk's reference temperature of 100 °C.

**What fixes the numbers.** Zhou 2024 fed 20 mmol/L alanine with 20 mmol/L glyoxal or methylglyoxal
in water at initial pH 8 and 100, 110, 120 °C and printed the pyrazine and 2,5-dimethylpyrazine
formation rates (0.0279 / 0.0791 / 0.1507 and 0.0035 / 0.0100 / 0.0230 µmol L⁻¹ min⁻¹) with the
barriers 100.6 and 111.7 kJ/mol. The conversion is below 0.2 %, so each printed "zero-order" rate is
an initial rate at 20 mM × 20 mM and re-expresses as a second-order constant: at 100 °C
6.98 × 10⁻⁸ L/(mmol·min) for the glyoxal route and 8.75 × 10⁻⁹ for the methylglyoxal route
(`zhou2024_extraction.md` section 4). Those six rates are the FIT rows for `k_go_pz` and
`k_mgo_dmp` (log10 k at 100 °C and Ea each), with the printed barriers as the priors' centres and
the three-point refit discrepancy (100.6 to 103.1; 111.7 to 114.9) as their band.

**Two declared transfers, each stated on the card.** (i) Alanine to glycine: the ring carbons come
from the dicarbonyl and the amino acid supplies the nitrogen, so the products are the same; the rate
is not. No paper measures the pair on a fed dicarbonyl. Leahy 1989 shows the amino acid matters at
the whole-cascade level (lysine over asparagine 35× for pyrazine, 6× for methylpyrazine). The
transfer is declared with a ±0.5 dex band, the largest declared uncertainty of the wave. (ii) The
mixed route: no rate exists for glyoxal + methylglyoxal to methylpyrazine. `k_mixed_mp` is declared
as the geometric mean of the two measured constants with the two as its band, and methylpyrazine
is reported, not shipped, unless T3 passes for it.

**pH.** Zhou's pots are unbuffered at pH 8; the trunk's pots sit at pH 5 to 7. Leahy & Reineccius
1989 (chapter 18) give the same laboratory's lysine + glucose pyrazine rate at pH 9, 7 and 5
(3.596 / 1.346 / 0.0938 ppm/h at 95 °C; methylpyrazine likewise), i.e. k(9)/k(7) about 2 to 3 and
k(9)/k(5) about 20 to 60. Within-study ratios may be fitted. The step gets one pH term, a
piecewise-linear log10 k against pH with a free slope above 7 and a free slope below 7, shared by
the three routes, fitted on the four Leahy ratios (pyrazine and methylpyrazine, pH 7 and pH 5
against pH 9). The buffer changes between the arms (borate at 9, citrate-phosphate at 7 and 5), and
that confound is recorded, not corrected.

## 3. What runs

A generator derived from the trunk's shipped wave, with the three reactions and the six new
coordinates (`k_go_pz`, `k_mgo_dmp` as log10 k at 100 °C and Ea; two pH slopes) added the way B10
added the route barriers; every existing trunk coordinate frozen at its shipped value. Fit rows: the
six Zhou rates and the four Leahy pH ratios (ten rows, six coordinates; the methylpyrazine rows
inform only the pH slopes). Two starts, the 600-evaluation budget, then the Laplace covariance.

Hold-outs, never in the objective:
- Leahy 1989 chapter 7, lysine + glucose at 95 °C, pH 9, 2 h: total pyrazines 13.1 ppm and the
  Table III distribution (pyrazine : methylpyrazine : 2,5-dimethylpyrazine), with glycine standing
  in for lysine (declared) and the same 0.1 M + 0.1 M charge.
- Leahy 1989 chapter 7, the barriers 150 / 153 / 177 kJ/mol (whole cascade, 75 to 95 °C): a band the
  model's apparent barrier from a glucose + glycine pot is compared with, not a fit.
- Yu 2018, thermal arm: the barriers 99.8 ± 6.7 (2,5-dimethylpyrazine) and 104.1 ± 5.1
  (tetramethylpyrazine) kJ/mol on a glucose + glycine pot at pH 10, 70 to 90 °C; the rate constants
  themselves are not used (their time unit is not printed).
- Zhou 2023's pyrazine and methylpyrazine columns (cysteine + xylose Amadori, 120 °C, pH 6 to 8) as
  an out-of-lane report: the sulfur lane holds methylglyoxal and cysteine, and cysteine's Strecker
  route is not this wave's.
- The two `_Internal2026` isolate bundles that carry a 2,5-dimethylpyrazine value are synthetic
  legacy output and are not hold-outs.

## 4. What counts as success, declared before the run

- **T1, the fit rows.** Each of the six Zhou rates within 0.3 dex; each fitted barrier inside its
  printed-to-refit band.
- **T2, nothing else moves.** Every scored panel row changes by less than 0.05 dex; the pyrazine
  flux is 0.2 % of the dicarbonyl in Zhou's pot and must stay a spectator to the trunk.
- **T3, another laboratory's shape.** Leahy's 95 °C distribution pyrazine : 2,5-dimethylpyrazine
  within 0.5 dex, glycine for lysine declared. Methylpyrazine's share is reported; if it is within
  0.5 dex the mixed route ships, else methylpyrazine stays refused with the reason on the card.
- **T4, another laboratory's level.** Leahy's 2 h total within tenfold (validation; reported).
- **T5, identification.** Laplace sigma below one decade on the four rate coordinates; the two pH
  slopes reported with theirs (four ratios may not identify both; if not, the lower slope is held
  at Leahy's pH 9 to 5 ratio and the card says so).
- **T6, the pH direction.** The fitted k(pH 5)/k(pH 9) lies between 1/60 and 1/20.

Ships as B18 if T1, T2 and T5 pass; T3, T4 and T6 are reported. If T2 fails the step has a
bookkeeping error, not a chemistry result, and the wave is fixed and re-run under the same
pre-registration.

## 5. What it will not do

It will not make alanine-specific ethylpyrazines, tetramethylpyrazine (which needs the
acetoin / diacetyl route), methoxypyrazines, or the sulfur lane's cysteine-Strecker pyrazines; it
will not move any shipped trunk coordinate; it will not read a hold-out. The parent pyrazine needs a
molecule row in the compound registry before its prediction can be keyed (the registry today has
only the class alias); that is part of the wave's engineering, through the registry's generator.

## 6. Outcome (2026-09-08, run the same evening) — SHIPS, with two conditionalities on every answer

**What was built.** Five steps rather than three: the pre-registered net reactions were fourth order
in mass action (two dicarbonyls and two glycines), which is not Yu 2018's rate law, so the step was
written as chemistry says it goes: two Strecker deaminations (dicarbonyl + glycine → aminoketone +
CO₂ + formaldehyde), second order and rate-determining, and three aminoketone condensations on one
constant declared fast (Jousse 2002 call the condensation "fast"; Zhou's products rise linearly from
the start). With the condensation fast the aminoketones sit at steady state and the pyrazine rate is
the Strecker rate over two; the mixed pyrazine follows the two pools statistically, which is the
pre-registered "geometric mean" with its factor two. Methylpyrazine is keyed `MPZ` (MP is a sulfur
species). The pH term is stored at the trunk's reference pH 6.8, so every condition term is 1 at the
references; Zhou's pH-8 pots carry its factor in the fit. Generator
`generate_kinetic_core_b18_fit.py` (its own ten-row objective on the trunk integrator, two starts,
Laplace inside); ship rule `generate_kinetic_core_b18_ship_rule.py`.

**The fit.** Cost 2.64 on 10 rows, 6 free (reduced χ² 0.66); both starts agree to 2 × 10⁻⁸. The six
Zhou rates within 0.07 dex, the four Leahy ratios within 0.13 dex. log10 k at 100 °C, pH 6.8:
glyoxal route −6.54, methylglyoxal route −7.53 (Laplace σ 0.08 and 0.08 dex); barriers 103.1 and
114.9 kJ/mol, both on the upper edge of their printed-to-refit band (σ 16 kJ/mol: a three-point
ladder cannot pin a barrier inside a 3 kJ/mol band); pH slopes 0.197 above 7 and 0.580 below
(σ 0.04 and 0.05).

- **T1 passed**: worst Zhou row +0.071 dex; barriers inside their bands.
- **T2 passed**: 39 predicted panel numbers compared with the tracked scorecard; the largest change is
  5 × 10⁻¹⁶ dex. The pyrazine flux is a spectator to every scored row.
- **T3 failed, and says something.** Leahy's 95 °C / 2 h distribution is pyrazine : 2,5-dimethyl-
  pyrazine 23 : 1; the model gives 6 × 10⁻⁸ : 1, because the trunk makes glyoxal only through the
  B13 dry-glass entry (glucose → glucosone, the smallest constant on the trunk) and so makes almost
  no pyrazine from a sugar + amine pot. Methylpyrazine's share is off by 3.7 dex the other way; the
  mixed route is reported, not shipped as a prediction of its own.
- **T4 failed**: the 2 h total is 2.9 decades low (17 µg/L against 13.1 mg/L), glycine for lysine.
- **T5 passed** on the two log10 constants (σ 0.08 dex); the barriers sit on the band edge as said.
- **T6 passed**: k(pH 5)/k(pH 9) = 0.0257, inside 1/60 to 1/20.
- **Barriers against the hold-out laboratories**: the model's apparent barriers from a glucose +
  glycine pot are 390 kJ/mol (Yu 2018's 2,5-dimethylpyrazine at pH 10, 70 to 90 °C; measured
  99.8 ± 6.7) and 300 to 460 kJ/mol (Leahy at pH 9, 75 to 95 °C; measured 150 to 177). The step's own
  barriers are 103 and 115; the excess is the trunk's dicarbonyl supply, whose temperature dependence
  in water at 70 to 95 °C the model gets steeply wrong.

**The directional panel.** Two 2,5-dimethylpyrazine claims that the engine refused for want of the
species now evaluate, and both agree: PH-06 (the compound rising with pH from 4 to 9) and TEMP-04
(160 °C over 145 °C in an extruded pea-protein system). The headline moves from 23 of 41 to 25 of 43.

**Ships, by the rule (T1, T2, T5).** The frozen literals are in `parameters_pyrazine.py`; a unit test
asserts they equal this report. Two conditionalities travel on every pyrazine answer as warnings:

1. **The supply is not measured.** The two Strecker constants are measured on fed dicarbonyls; from a
   sugar + amine pot the yield follows the trunk's dicarbonyl supply, which T3 and T4 show to be
   wrong by orders of magnitude at 95 °C in water. A pyrazine number from this model is a
   fed-dicarbonyl statement.
2. **The glyoxal sink.** The B13 dry-glass glyoxal sink (180 °C, barrier fixed to zero) removes 98 %
   of a fed 20 mM glyoxal in two hours at 100 °C, so the modelled Zhou pot's pyrazine growth is not
   linear (0 to 60 min rate 1.75 × the 0 to 120 min rate) where Zhou's is; with that sink zeroed
   (`--variant nosink`, information only) the glyoxal Strecker constant is 0.59 dex lower and the
   growth linear. The methylglyoxal pot loses its dicarbonyl too (93 to 98 %), through the trunk's
   fitted melanoidin sink and the measured furanone step, which this wave may not move.

**What this asks for next.** A wave on the small dicarbonyls in water: their formation from a sugar
+ amine pot at 70 to 120 °C (Leitzen 2021's aqueous glucose series is on disk and already misses,
DIC-01), and their loss (Zhou 2024's glyoxal and methylglyoxal time courses are in Figure 4,
figure-only; the wishlist's "glyoxal loss at two temperatures" stands). Until then the pyrazine
answer's first line is its own caveat.

**Post-run note (2026-09-08, later the same night).** Zhou et al. 2025 (`zhou2025b_extraction.md`),
the same laboratory, prints alanine + glyoxal pyrazine formation rates at 70 / 80 / 90 °C
(0.1165 / 0.2806 / 0.8867 µmol L⁻¹ min⁻¹, Ea 105.0 kJ/mol) whose Arrhenius line has the 2024 slope
within 4 kJ/mol and an intercept seventy times higher: the 2024 line predicts 0.012 at 90 °C where
2025 prints 0.887, and the 2025 paper never states the reactant concentrations of those runs. The
two ladders cannot both be initial rates at 20 + 20 mM. B18 was fitted on the 2024 ladder, whose
concentrations are printed; the discrepancy is recorded here and in the backlog, and until the
laboratory's concentrations are known it is a third conditionality on the glyoxal route (up to
1.85 dex). The 2025 paper also gives the first rate for alanine + xylose → the Amadori compound
(0.0034 / 0.0060 / 0.0170 mmol L⁻¹ min⁻¹ at 70 / 80 / 90 °C, Ea 83.1 kJ/mol), a trunk quantity the
sulfur lane's pentose Amadori step could be checked against later.

