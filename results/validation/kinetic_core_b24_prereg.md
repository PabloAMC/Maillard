# Pre-registration: wave B24, 2-acetyl-1-pyrroline from proline (written 2026-09-09, before the fit ran)

*(B23 is not used as a wave name: `generate_kinetic_core_b2_3_fit.py` is imported as `B23` throughout
the sulfur generators.)*

## 1. Why

2-Acetyl-1-pyrroline is the bread-crust and popcorn note of extruded and baked products, the one
proline odorant on the roadmap's coverage table, and no lane names it; until today it was not even in
the compound registry. Hofmann & Schieberle 1998b (`hofmann1998b_extraction.md`) fed 1-pyrroline with
methylglyoxal and printed yields at 100 °C in phosphate at pH 7 (fed-intermediate yields, which fit
under the owner's rule), fed proline with methylglyoxal at three ratios and printed the whole-chain
yields, and drew the mechanism: hydrated methylglyoxal acylates C-2 of 1-pyrroline, the aldehyde
carbon leaves as carbon dioxide, air oxidises the pyrrolidine to the pyrroline. Chan & Reineccius
1994 (`chan1994b_extraction.md`) give the only apparent barrier for the whole cascade in water, 60.2
kJ/mol. The trunk now carries methylglyoxal from Martins' aqueous step and the Strecker machinery of
B18; what is missing is proline and its two steps.

## 2. The arm

Three trunk-only species appended last: PRO, L-proline (C5 N1); PYRL, 1-pyrroline (C4 N1); AP,
2-acetyl-1-pyrroline (C6 N1). Two steps:

    r_mgo_pro   MGO + PRO  -> PYRL + 4 FRAG_C     k_mgo_pro   second order (the Strecker of proline: the ring nitrogen stays in 1-pyrroline; the dicarbonyl leaves as hydroxyacetone, booked to the fragment pool, and CO2)
    r_pyrl_ap   PYRL + MGO -> AP + FRAG_C         k_pyrl_ap   second order (the acylation; the aldehyde carbon of methylglyoxal as CO2)

Not written: the 2-acetyltetrahydropyridine branch (it needs hydroxyacetone as a species, and the
trunk has none), any loss of 1-pyrroline or of the product (no rate exists; Hofmann's excess-pyrroline
experiment shows a loss the arm cannot reproduce, and the ship rule reports that), and proline on
glyoxal (Hofmann fed methylglyoxal only). Proline is charged as PRO for these steps and, declared, as
glycine at the same molarity for the Amadori chemistry (the B22 declaration; proline is a secondary
amine and a weaker Amadori amine, so this is an upper bound on the dicarbonyl supply it makes).

**What is fitted.** Two coordinates, log10 `k_pyrl_ap` and log10 `k_mgo_pro` at 100 °C. Rows:
Hofmann 1998b Table 7 experiments 1 and 2 (1-pyrroline 2 mmol/L with methylglyoxal 10 and 2 mmol/L,
pH 7, 0.5 M phosphate, 100 °C, 30 min: 28.7 and 5.3 mol % of the pyrroline as the product) and Table
9 (proline 400 mmol/L with methylglyoxal 4, 40 and 400 mmol/L, same conditions: 0.0058, 0.0125 and
0.0179 mol % of the proline as the product). Five rows, each integrated as its pot with the fed
species charged directly and the yield read at 30 min; sigma 0.15 in log10 for the fed-pyrroline
rows, 0.25 for the proline rows (a whole chain with the Strecker supply inside it).

**What is declared.** (i) The barriers: none is measured for either step; the acylation carries
Chan & Reineccius 1994's whole-cascade 60.2 kJ/mol (the only 2-acetyl-1-pyrroline barrier on disk,
flagged as apparent), the proline Strecker carries B18's methylglyoxal Strecker barrier (114.9,
measured on glycine). (ii) Proline as glycine for the Amadori chemistry. (iii) No loss steps. (iv) The
pH term of B18 applies to the proline Strecker step (an amine-dependent step at pH 7 against the
trunk's 6.8 reference).

## 3. What runs

The generator integrates the five pots and compares the 30-minute yields with the printed ones; then
Hofmann's experiment 3 (1-pyrroline 10 + methylglyoxal 2 mmol/L: the model's yield per methylglyoxal
against the printed 0.33 mol %, where the source finds excess pyrroline SUPPRESSES the product); the
apparent barrier of the whole cascade in a glucose 100 + proline 100 mmol/L pot at pH 7 between 75
and 115 °C against Chan's 60.2 kJ/mol; and the B18 test pot's pyrazines before and after (proline
competes with glycine for methylglyoxal only when charged, so no change is expected).

## 4. What counts as success, declared before the run

- **T1, the rows.** The two fed-pyrroline rows within 0.3 dex; the three proline rows within 0.5 dex
  and in the printed order (the yield per proline rising with methylglyoxal).
- **T2, the panel is untouched.** No currently scored panel row moves by more than 0.05 dex.
- **T3, another laboratory.** The apparent barrier in the glucose + proline pot against Chan's 60.2
  kJ/mol: reported in kJ/mol (the arm's barriers are declared, so agreement would be luck).
- **T4, the suppression.** Hofmann's experiment 3: the model's yield per methylglyoxal against the
  printed 0.33 mol %; reported, expected to fail (no pyrroline loss is written).
- **T5, identification.** Laplace sigma below one decade on both coordinates, neither on its bound.

Ship rule: SHIP if T1, T2 and T5 hold; T3 and T4 are reported. If it ships, the answer for
2-acetyl-1-pyrroline on a proline pot carries the declared barriers and the no-loss note, and no
odour-activity value is computed for it (its threshold is not on disk). If it does not ship, the
steps stay at zero and the target is refused by name.

## 6. Outcome (2026-09-09, run the same day) — DO NOT SHIP, and the missing piece is named

**What was built.** The three species and two steps of section 2, trunk-only and appended last; proline
as PRO and, declared, as glycine for the Amadori chemistry; B18's pH term on the proline Strecker step;
2-acetyl-1-pyrroline added to the compound registry (threshold left null, not on disk); generator
`generate_kinetic_core_b24_fit.py`, ship rule `generate_kinetic_core_b24_ship_rule.py`.

**What happened.** The two fed-pyrroline rows fit: +0.18 and +0.30 dex (the acylation constant
2.74e-03 L/(mmol min) at 100 °C, three times Hofmann's own bilinear lower bound). The three proline rows do not:
-1.29, +0.13 and +1.16 dex at 4, 40 and 400 mmol/L of methylglyoxal. The model's whole-chain yield rises
almost linearly with the methylglyoxal charge (a thousandfold over the ladder) where Hofmann's rises
threefold, because in the source the excess amino acid converts the methylglyoxal to hydroxyacetone
and 1-pyrroline that then go to the tetrahydropyridine (the switch the authors name, Table 9's
AP : ATHP ratio 0.16 to 12.8), and the arm has neither that branch nor any loss of 1-pyrroline: in
experiment 3 (1-pyrroline in fivefold excess) the model makes 41 mol % of the methylglyoxal into the product
against the printed 0.33 (+2.1 decades). Both coordinates are identified (sigma 0.5 and 0.7 decades), the panel does not
move, and the apparent barrier of the cascade in a glucose and proline pot is 591 kJ/mol against Chan's 60 (the
trunk's dicarbonyl supply, not the step). T1 fails; verdict by the rule: DO NOT SHIP.

**What it teaches.** The acylation step is measured well enough to carry; the chain from proline is
not a chain of two steps. A next pre-registration (B24b) needs hydroxyacetone as a trunk species (the
Strecker of an amino acid on methylglyoxal makes it; Hofmann's Table 4 gives its reaction with
1-pyrroline at four pH values, Table 5 the intermediate's own conversion), the tetrahydropyridine as
the competing product, and a loss of 1-pyrroline sized on experiment 3. With those, Table 9's switch
is the within-study shape to fit. The steps stay in the network at zero, the optimum is recorded in
`parameters_proline.FROZEN_B24`, and a request for 2-acetyl-1-pyrroline is refused with this verdict.

