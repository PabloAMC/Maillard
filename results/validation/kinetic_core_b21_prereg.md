# Pre-registration: wave B21, the aqueous glucosone route to glyoxal (written 2026-09-09, before the fit ran)

## 1. Why

The pyrazine step (B18) ships with the caveat that from a sugar and amino acid pot the model makes
far too little glyoxal in water: the total pyrazine misses Leahy 1989's 95 °C pot by 2.9 decades,
and the miss is the supply, not the step (`kinetic_core_b18_prereg.md`, section 6). The trunk's
only route to glyoxal is glucose → glucosone → glyoxal with constants from an amine-free sugar
glass at 160 to 200 °C (B13, Kocadagli 2016): at 100 °C the entry runs at a millionth per minute.
Four dossiers read on 2026-09-08 and 2026-09-09 (`hamzalioglu2026`, `xia2022`, `quan2020`,
`yu2020`) say what happens in water: the glucosone comes from the Amadori compound, not from the
sugar (Hamzalioglu 2026 fits lactulosyl-lysine → glucosone at 110 to 140 °C with a barrier, the
only aqueous entry on disk); glyoxal is made "far more" than methylglyoxal at 130 °C (Xia 2022);
and a glucose and lysine pot holds about 0.1 mmol/L of glyoxal at 100 °C and 0.6 at 130 °C within
twenty minutes (Quan 2020, levels).

## 2. The step

One new step on the trunk lane, on existing species:

    r_ama_g    AMA -> G + Gly        k_ama_g    first order (the Amadori compound's oxidative cleavage to glucosone, returning the amine)

and the glucosone → glyoxal step already on the trunk (`r_g_go`, `k_g_go`) given an aqueous value.
Carbon balances: the Amadori compound (C8 N1) gives glucosone (C6) and glycine (C2 N1). The dry-glass
glucose → glucosone entry (`k_glc_g`) and the dry-glass glyoxal sink (`k_go_sink`) are left as they
are: the sink's effect is B18's declared conditionality (0.59 decades on the glyoxal Strecker
constant) and changing it would move a shipped fit; its removal is reported as a variant, not
shipped.

**What fixes the numbers.** Hamzalioglu 2026 Table 1 (whole milk, lactose and casein-bound lysine,
110 / 120 / 130 / 140 °C, multiresponse fit): lactulosyl-lysine → glucosone 1.9e-2, 2.1e-2, 2.5e-2
(± 7.1e-2), 1.5e-1 per minute, barrier 75.9 ± 21.1 kJ/mol; glucosone → glyoxal 7.4e-2 (± 1.5),
3.3e-1, 3.5e-1 (± 3.7e-1), 2.4e-3 (indeterminate) per minute, barrier 4.2 ± 15.7 kJ/mol. First-order
constants need no water basis, so they transfer as printed. Fit rows: the four glucosone-formation
constants (the 130 °C one weighted by its wide interval) and the two determinate glyoxal-formation
constants (120 and 130 °C); the other two reported. Six rows.

**What is fitted.** Two coordinates: log10 `k_ama_g` at 100 °C and log10 `k_g_go` at 100 °C, bands
two decades either side of the prior centre (Hamzalioglu's values brought to 100 °C with the
declared barriers).

**What is declared.** (i) The barriers: Hamzalioglu's measured 75.9 kJ/mol for the glucosone
formation; for glucosone → glyoxal the measured 4.2 kJ/mol is consistent with zero and the trunk's
glass value is 93.8, so the aqueous constant carries Hamzalioglu's barrier and the glass barrier
is reported as the other reading. (ii) The transfer from lactulosyl-lysine in milk to
fructosyl-glycine in water: the oxidative cleavage is on the sugar moiety and the amine leaves, so
the products are the same; the rate is not measured for the pair; a ± 0.5 dex band is declared on
every glyoxal, glucosone and pyrazine answer (the B18 precedent). (iii) The lysine of Quan 2020's pot
stands in as glycine at the same molarity (declared, as Leahy's lysine did in B18).

## 3. What runs

The generator compares each fitted constant, evaluated at the row's temperature with its declared
barrier, with the printed value (no integration in the objective), then integrates three pots:
Quan 2020's (glucose 100 + glycine-for-lysine 30 mmol/L, 0.1 M phosphate pH 7, 100 and 130 °C, 21
min) for the glyoxal level; Xia 2022's (200 + 200, pH 7.5, 130 °C, 80 min) for the glyoxal to
methylglyoxal ordering; Leahy 1989's (100 + 100, pH 9 borate, 95 °C, 2 h) for the total pyrazine
B18 missed; and re-scores the B1 browning hold-out (Martins 2005's melanoidin trajectory at 80, 100
and 120 °C) with the new step in the network, through the frozen B1 hold-out generator's own
scoring function.

## 4. What counts as success, declared before the run

- **T1, the rows.** The four determinate rows (glucosone formation at 110, 120 and 140 °C; glyoxal
  formation at 120 °C) within 0.3 dex; the wide ones reported.
- **T2, nothing shipped breaks.** (a) The B1 browning hold-out's median fold error stays below 2
  (shipped 1.45) and its fraction within threefold stays 1.0; (b) no currently scored panel row
  moves by more than 0.3 dex; the exact changes are listed.
- **T3, the glyoxal level.** Quan 2020's pot: the model's glyoxal at 21 min within the printed range
  widened by 0.5 decades either side, at 100 °C (0.052 to 0.127 mmol/L) and at 130 °C (0.144 to
  0.605).
- **T4, the ordering.** Xia 2022's pot at 130 °C: glyoxal above methylglyoxal at 80 min.
- **T5, the pyrazine supply.** Leahy's total pyrazine at 95 °C: the model's miss, in decades, against
  B18's 2.9; reported (the step is glycine's, Leahy's is lysine).
- **T6, identification.** Laplace sigma below one decade on both coordinates, neither on its bound.

Ship rule: SHIP if T1, T2, T3 and T6 hold; T4 and T5 are reported. If it ships, the frozen literals
enter `parameters_dicarbonyl.py` beside the glass constants, every glyoxal and pyrazine answer
carries the transfer band, and B18's supply caveat is rewritten to what remains. If T2a fails, the
step is DO NOT SHIP however well the rows fit: the browning hold-out is the trunk's one
out-of-sample success and the step must not spend it.

## 6. Outcome (2026-09-09, run the same day) — SHIP

**What was built.** The step `r_ama_g` on the trunk and the aqueous value of `k_g_go` in the
operative set, the glass value kept in `parameters_dicarbonyl.DICARBONYL_PARAMETERS` as the record;
`FROZEN_B21` asserted equal to the report by `tests/unit/test_kinetic_core_b21.py`; the transfer
caveat on every glyoxal, glucosone and pyrazine answer. Generator `generate_kinetic_core_b21_fit.py`,
ship rule `generate_kinetic_core_b21_ship_rule.py`.

**What happened.** Cost 3.54 on six rows (reduced chi-square 0.89); both coordinates
identified (sigma 0.09 and 0.14 decades), neither on a bound. T1 passes: the decisive rows within
0.18 dex; the wide 130 °C glucosone row sits 0.36 dex above its printed centre, inside its interval.
T2 passes on both counts: the B1 browning hold-out's median fold error moves from 1.43 to 1.31 with
every point still within threefold (the new route takes a share of the Amadori flux that the browning
did not need), and not one scored panel row moves (the panel's trunk rows are HMF and browning, which
the glucosone route does not touch at the 1e-9 level). T3 passes: in Quan 2020's pot the glyoxal at 21
min is 0.023 mmol/L at 100 °C against the printed 0.052 to 0.127 (-0.35 dex from the range) and 0.387 at
130 °C, inside 0.144 to 0.605; before the wave the same pot held 3.2e-06 and 4.3e-04 mmol/L, four and two and a
half decades low. T4, reported: at 130 °C and 80 min in Xia's pot the model's glyoxal is 1.23 mmol/L and
its methylglyoxal 7.82, so the ordering Xia reports (glyoxal far above) is not reproduced: the model's
methylglyoxal from Martins' aqueous step is high, or Xia's derivatisation over-reads glyoxal (their
own flag). T5, reported: Leahy's total pyrazine moves from -2.88 to -2.78 decades. Verdict by the rule: SHIP.

**What it changes, and what it does not.** A sugar and amine pot now holds glyoxal at the level a
laboratory measures, and pyrazine (the parent) is no longer absent: in the B18 test pot the order of
the three pyrazines turns from 2,5-dimethyl > methyl > parent to parent > methyl > 2,5-dimethyl, which
is Xia's direction for the dicarbonyls. The pyrazine step's supply caveat is rewritten accordingly.
What it does not change: Leahy's total. The 2.8 decades that remain are not glyoxal; they sit in the
Strecker step at 95 °C and pH 9, in the lysine that stands in as glycine, or in Leahy's own recovery,
and the caveat now says so. The B13 glyoxal sink is untouched, so B18's conditionality stands.

**An unforeseen finding, recorded after the verdict.** The pre-registration checked the browning
hold-out and the panel; it did not check B1's own fit rows, Martins 2005's nine concentration
series at 80, 100 and 120 °C. Re-scored with the frozen B1 vector and the live network, the half
sum of squares rises from 2017 to 2495 (+24 %); the Amadori compound's median error goes from 0.035
to 0.093 dex (the model now holds a fifth less of it at 100 °C and two hours), formic acid from 0.089
to 0.123, while 3-deoxyglucosone, methylglyoxal and 1-deoxyglucosone move by less than 0.01 dex. The
new route drains the Amadori compound at nine thousandths per minute at 100 °C, the size of Martins'
own Amadori → 3-deoxyglucosone step, and Martins' Amadori series does not want that drain at his
pH 6.8 and air-limited headspace. The verdict stands by the rule as written; the tension is stated
on every glyoxal answer's caveat, the B1 regression pin is moved with the arithmetic written out,
and the next pre-registration on this route (B21b) is a joint fit of the glucosone rate against
Hamzalioglu's constants AND Martins' Amadori series, which will either find a rate both accept or
show that the milk constant does not transfer.

