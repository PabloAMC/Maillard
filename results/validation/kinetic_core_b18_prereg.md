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

## 6. Outcome

*(not yet run; generator to be derived from the trunk's shipped wave as described in section 3)*
