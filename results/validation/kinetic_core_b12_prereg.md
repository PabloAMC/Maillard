# Wave B12 pre-registration — water activity and pH on the trunk lane (programme step W4, part 1)

*Written 2026-09-07 before the directional panel was re-scored. Owner's instruction: "the trunk lane
from the eleven unread kinetics papers ... would gain water activity and the dicarbonyls." B12 is the
first of two trunk waves: the two condition terms. B13 (dicarbonyls: glucosone, glyoxal, diacetyl;
the Lee 2022 / 2024 and Brands 2002 validation sets) follows. Module: `src/kinetic_core/trunk_conditions.py`;
licence: FIT_HOLDOUT_DECLARATION.md Amendment 21.*

## 1. What the eleven papers turned out to hold, for the trunk

| paper | what it gives the trunk | used in B12 |
|---|---|---|
| Pereyra Gonzales 2010 (milk powder, lysine loss, a_w 0.33-0.98, 37-60 C, 95 % CI) | the NET a_w shape of the Maillard rate at fixed dry-basis composition: ~3.5x at a_w 0.4-0.7 over solution, falling toward the glass | **the a_w multiplier's centre** |
| Bell 1995 (PVP / glucose / glycine, fixed molality, a_w 0.33-0.96, 25 C) | deconfounded: flat above the glass, ~5x lower in it; reactant molality dominates | **the band's floor (no effect)** and the glass warning |
| Lievonen 2002, Miao 2004 | browning in amorphous carbohydrate glasses; Ea rises as water falls (156 -> 118 kJ/mol) | recorded: the a_w effect is temperature-dependent and B12 does not model that |
| Martins & van Boekel 2003 Part II (DFG degradation, 100/120 C x pH 5.5/6.8, HPD) | every Amadori-decay step 3-16x faster at pH 6.8: 0.37-0.92 decades per pH unit | **the pH exponent (0.69) and its band** |
| Kocadağlı & Gökmen 2016 JAFC (glucose glass, 160-200 C) | measured rates + barriers for the amine-free caramelisation steps | ALREADY in the trunk since B7 (the census had this wrong: only the Food Chem wheat-flour paper is unread) |
| Kocadağlı 2016 Food Chem (glucose / wheat flour), Sen 2022 (nuts), Ağçam 2022 (juices), Gürsul Aktağ 2020 (juice dicarbonyls) | low-moisture and acid-route validation sets; glucosone / glyoxal / diacetyl kinetics | B13 |
| Hidalgo 1993, Zamora 2013 | lipid-carbonyl browning and Strecker barriers | lipid lane / priors; not trunk |
| Göncüoğlu Taş 2017 | structure only, no barrier | none |

## 2. What changes

Two DECLARED terms, no fit: measured within-study ratios installed as constants with bands, the same
standing as a measured barrier override.

- **a_w multiplier on the amine-sugar condensation** (`k_schiff`, and through the pinned split the
  Amadori rearrangement): piecewise-linear through Pereyra Gonzales's k(a_w)/k(0.98) at 50 and 60 C
  (3.05 / 3.73 / 3.40 / 3.68 / 2.58 / 1.00 at a_w 0.33 / 0.43 / 0.52 / 0.69 / 0.85 / 0.98); band from 1.0
  (Bell's fixed-molality plateau) to 1.2x the centre; held at the 0.33 value below the table with a
  glass warning. Dehydration steps carry no a_w term (declared gap).
- **pH factor on the three Amadori-decay steps** (`k_ama_tdg`, `k_ama_odg`, `k_ama_mgo`):
  10^(0.69 (pH - 6.8)), band exponent (0.37, 0.92), measured window 5.5-6.8, extrapolation warning
  outside it. The amine-free caramelisation entries carry no pH term (no source contrasts them).
- **Engine:** the trunk lane applies both before integration (`_integrate_program`); the envelope
  declaration prints them; `axis_refusal` answers a_w moves when a resolved lane carries the term
  (trunk) and pH moves on the trunk. The acrylamide lane keeps both refusals (De Vleeschouwer 2009,
  on disk and unread, is its a_w source; a B13 item). The copy of the trunk steps inside the sulfur
  network is untouched (a sulfur wave's decision).
- **At the references (a_w None or >= 0.98, pH 6.8) every factor is exactly 1.0**, so B1, B7 and every
  panel row are reproduced bit for bit (the trunk's one panel row, Steinhagen 2021, is amine-free at
  a_w 0.99: no Amadori flux, no change).
- **The bands are printed, not yet sampled:** the envelope's draw of the multiplier and the exponent
  is a B13 item (`CoreDraw` hooks exist in `trunk_conditions.apply`).

## 3. Pre-registered expectations on the directional panel

`AW-01` (HMF falls with water content, a_w 0.3 / 0.6 / 0.9, glucose + glycine, 140 C): with a
condensation multiplier of 3.05 / 3.6 / ~1.7 the model will predict HMF LOWER at 0.3 than at 0.6 --
**expected DISAGREE**, and honestly so: Pereyra Gonzales's shape has a maximum, not a monotone fall.
`AW-03` (Maillard rate peaks at a_w 0.6-0.7, falls at both 0.25 and 0.9): **expected AGREE** if the
0.9 arm falls below the 0.65 arm by more than the 5 % flat tolerance (multiplier ~1.7 vs 3.7).
`AW-02` (acrylamide, bell-shaped): stays not evaluable (acrylamide lane). pH claims on the trunk:
`CYS-02` and any HMF pH claim become evaluable; no expectation is pre-registered for them (the
term acts on Amadori decay only, and HMF in those pots comes through both routes).

## 4. What would falsify the terms

A within-study a_w series on a trunk observable (HMF, 3-DG, browning) at fixed dry-basis
composition that does NOT rise from solution toward a_w 0.5-0.7 by at least 2x at 40-120 C; or a
pH series on the Amadori branch with the opposite sign. Either is a measurement, not a refit.

## 5. Outcome (2026-09-07, after the re-score)

Both pre-registered expectations held. `AW-01` DISAGREE: HMF 1350 / 1500 / 629 ug/L at a_w 0.3 / 0.6 / 0.9,
a maximum rather than a monotone fall, as the Pereyra Gonzales shape implies (the claim's source, Ma 2024,
is an extrusion study whose "moisture" is also a temperature and residence-time change). `AW-03` AGREE:
1350 / 1530 / 762 ug/L at 0.25 / 0.65 / 0.9, a peak. `AW-02` (acrylamide) stays not evaluable. `CYS-02`
(cysteine present vs absent, HMF) is now answered on the trunk and still disagrees (629 vs 674). Panel
headline 17/28 -> **18/30**; pH-and-water-activity 4/5 -> 5/7; independent not-evaluable 29 -> 27; the
absolute scorecard is unchanged (4/39, 3/38). The engine's refusal list for water activity now names the
lane, and the wishlist's "blocked: no lane carries a water-activity term" line is gone for the trunk.
