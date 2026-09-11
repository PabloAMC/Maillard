# Pre-registration: wave B22, the methionine chain on the sugar path (written 2026-09-09, before the fit ran)

## 1. Why

Methional is the cooked-potato and meaty note that sits at the top of the desirable list, and
neither it nor its children (methanethiol, dimethyl disulfide, dimethyl trisulfide) is named by
any lane: the coverage figure shows "nothing" for the first three. The corpus now holds the first
rate table for the chain (`pan2025_extraction.md`: methional, methanethiol, dimethyl disulfide and
dimethyl trisulfide accruing at zero order in a methionine and fruit-sugar pot at 100, 120 and
140 °C), a level series from a second laboratory (`deng2022_extraction.md`: methional against time
at 120 °C from methionine and glucose, and from the methionine Amadori compound), and the
oxidant-resolved fate of methanethiol from a third (`chin1994_extraction.md`). The Strecker step
that makes methional is the one the pyrazine step (B18) already carries for glycine on glyoxal and
methylglyoxal, and B21 has just given the trunk an aqueous glyoxal supply; what is missing is the
amino acid's identity and the two steps after the aldehyde.

## 2. The arm

Four trunk-only species appended last: MET, methionine (C5 N1 S1); MTAL, methional (C4 S1); MSH,
methanethiol made from methional (C1 S1; the sulfur lane's own methanethiol, MESH, from thiamine,
is a different pool on a different lane and keeps its name); DMDS, dimethyl disulfide (C2 S2).
Four steps:

    r_go_met     GO  + MET -> AKG + MTAL + FRAG_C      k_go_met     second order (Strecker on glyoxal; the aminoketone is glycine's, the aldehyde is methionine's)
    r_mgo_met    MGO + MET -> AKM + MTAL + FRAG_C      k_mgo_met    second order (Strecker on methylglyoxal)
    r_mtal_msh   MTAL -> MSH + 3 FRAG_C                k_mtal_msh   first order (the retro-Michael release of methanethiol; acrolein to the fragment pool)
    r_msh_dmds   2 MSH -> DMDS                         k_msh_dmds   second order, APPARENT: the pot's internal oxidant is not tracked on this lane

Dimethyl trisulfide is not written: it needs hydrogen sulfide (Chin & Lindsay 1994 detect none
without it), which the sugar path does not carry; a request for it is refused with that reason.
Methionine is charged as MET for the Strecker step and, declared, as glycine at the same molarity
for the Amadori chemistry that makes the dicarbonyls (the same amine plays both roles in turn; at
Pan's 0.27 mmol/L the Amadori consumption of amine is negligible within ten minutes). The
aminoketones the Strecker step makes are glycine's AKG and AKM by construction (the ring carbons of
a pyrazine come from the dicarbonyl), so methionine feeds the pyrazine pools as glycine does.

**What is fitted.** Four coordinates: the identity ratio of methionine to glycine on the two
Strecker steps, one log10 factor applied to both B18 constants (`k_go_met = r × k_go_ak`,
`k_mgo_met = r × k_mgo_ak`, B18's barriers and pH term kept); log10 `k_mtal_msh` at 100 °C and its
barrier within 20 to 150 kJ/mol; log10 `k_msh_dmds` at 100 °C with its barrier declared from Pan's
apparent value for dimethyl disulfide (81 kJ/mol; the pair of methanethiol and disulfide rates is
what identifies it). Rows: Pan 2025 Table 2's nine zero-order constants (methional, methanethiol,
dimethyl disulfide at 100, 120, 140 °C), read as micromoles per litre per second (the unit the
paper does not print; the two printed endpoint levels fix it, dossier flag 1), each modelled as the
mean formation rate over 30 to 600 s in Pan's pot (methionine 0.268 mmol/L, fructose 111 + glucose
83 mmol/L, sucrose 44 omitted and declared, 50 mmol/L citrate pH 6.2), sigma 0.15 in log10 for
methional and methanethiol, 0.25 for the disulfide (its regressions are weaker). Two starts, scipy
least squares, Laplace at the optimum.

**What is declared.** (i) The unit of Pan's constants (inferred; the wave stands or falls with it,
and says so). (ii) Sucrose omitted: the trunk has no sucrose, and Pan's sugars are 4 % consumed in
four minutes at 120 °C. (iii) The identity ratio is conditional on the trunk's dicarbonyl supply in
Pan's pot (B21's glyoxal, Martins' methylglyoxal): a fed-dicarbonyl methionine experiment would free
it, and the wishlist asks for one. (iv) The disulfide step's constant is apparent and pot-specific
(Xu 2010: the oxidant is internal); it carries a flag and no claim of transfer. (v) Methionine as
glycine for the Amadori chemistry.

## 3. What runs

The generator integrates Pan's pot at the three temperatures with the candidate constants and
compares the mean formation rates with the printed ones; then Deng 2022's pot (methionine 200 +
glucose 200 mmol/L, initial pH 7.5 unbuffered, 120 °C) for the methional level at 30, 60, 120 and
180 min against the printed 18.5, 25.3, 75.8, 87.8 µg/L (response factor 1, so an order of
magnitude); and the B18 test pot for the three pyrazines before and after (methionine feeds the
aminoketone pools).

## 4. What counts as success, declared before the run

- **T1, the rows.** The six methional and methanethiol rows within 0.3 dex; the three disulfide
  rows within 0.5; the fitted barrier of the release step off its bounds.
- **T2, nothing shipped breaks.** No currently scored panel row moves by more than 0.05 dex (no
  panel pot charges methionine, so the arm is inert there).
- **T3, another laboratory.** Deng 2022's methional at 120 min within one decade (the source's
  response factor is 1, so tighter is not claimed) and rising from 30 to 120 min as the source's does.
- **T4, identification.** Laplace sigma below one decade on every coordinate; the identity ratio
  within its band of ± 2 decades around 1.
- **T5, the fate.** Chin & Lindsay 1994's methanethiol half-life with copper at 30 °C (about 17
  minutes) against the model's disulfide step extrapolated to 30 °C at their methanethiol level:
  reported, in decades, with the note that theirs is an added oxidant.

Ship rule: SHIP if T1, T2 and T4 hold; T3 and T5 are reported. If it ships, `explain methional`
answers with the arm and its three declarations, the coverage figure moves methional, methanethiol
and dimethyl disulfide to "modelled", and every methional answer carries the unit and the supply
conditionality. If it does not ship, the steps stay at their frozen literals with a DO-NOT-SHIP
note and the four targets are refused by name.

## 6. Outcome (2026-09-09, run the same day) — DO NOT SHIP, and the structure is refuted

**What was built.** The four species and four steps of section 2, trunk-only and appended last;
methionine charged as MET and, declared, as glycine for the Amadori chemistry; the pH term of the
pyrazine step applied to the two methionine Strecker steps; the fit generator
`generate_kinetic_core_b22_fit.py` (Pan's pot integrated at the three temperatures, the mean
formation rate over 30 to 600 s as the observable) and the ship rule
`generate_kinetic_core_b22_ship_rule.py`.

**What happened.** The optimiser drove the identity ratio to its ceiling, a hundred times glycine's
constants, and the disulfide constant to its ceiling too, and the methional rows were still
3.6 to 5.6 decades below Pan's printed rates
(methanethiol 2.0 to 5.3, the disulfide 2.6 to 10.6); cost 7374 on nine rows.
In the same run Deng 2022's pot (methionine 200 + glucose 200 mmol/L, 120 °C) came out
+1.5 to +2.9 decades ABOVE the printed methional and falling where the source rises. T1 fails, T4 fails
(three of four coordinates on a bound, the barrier of the release step at its floor), T2 passes (no
panel pot charges methionine). Verdict by the rule: DO NOT SHIP.

**Why, in the numbers.** In Pan's pot the trunk holds 9 µmol/L of glyoxal and 42 µmol/L of
methylglyoxal at 140 °C and ten minutes, against 268 µmol/L of methionine; with glycine's Strecker
constants (a few times 1e-7 litres per millimole per minute at 100 °C, from Zhou 2024's fed
dicarbonyls) the product of the three is a rate of nanomoles per litre per minute, and Pan measures
half a micromole per litre per minute. A hundredfold ratio does not bridge that; ten-thousandfold
would, and then Deng's pot, where the amine is a thousand times higher and the trunk makes tens of
milligrams of glyoxal per litre, overshoots by four decades instead of two. No single ratio serves
the two pots: methional does not form as free dicarbonyl times methionine with the Strecker
constants measured on fed dicarbonyls. Deng's own experiment says where it does form: the
methionine Amadori compound alone gives 1.4 to 2.6 times more methional than methionine plus
glucose, so the route is the Amadori compound's own decomposition (the sugar moiety supplies the
dicarbonyl in the same molecule), which is a first-order step in a methionine Amadori compound the
trunk does not carry. That is the structure the next pre-registration should write: MET + Glc → a
methionine Amadori compound (the trunk's Amadori formation with methionine as the amine), then
Amadori-Met → methional + the rest, with Deng's two series (methionine + glucose, and the fed
Amadori compound) as the rows, and Pan's rates as the second laboratory.

**What is kept.** The species and steps stay in the network at zero (the B17 precedent), the
optimum is recorded in `parameters_methionine.FROZEN_B22` and asserted against the report, and a
request for methional, methanethiol from methional or dimethyl disulfide is refused with this
verdict. Methionine as a precursor now resolves to the sugar path and is charged as glycine for the
Amadori chemistry, declared on every such answer; before this wave it was refused as unmapped.
Dimethyl trisulfide is refused with its own reason (no hydrogen sulfide on this path).

