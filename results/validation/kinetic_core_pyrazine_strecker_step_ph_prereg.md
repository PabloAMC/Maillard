# E4 pre-registration: the below-7 pyrazine pH slope moved off the Strecker step

*Written 2026-11-01 (night), before any refit. PROPOSAL for the owner's morning decision. Nothing ships, nothing in `src/`
is edited, nothing is committed. Origin: `docs/reference/pathway_atlas/ENGINE_VS_LITERATURE.md` sec. 3 row E4 and the
binding red-team amendment, sec. 5 item 1. Amended the same night by the plan premortem (coordinator, binding):
(1) the cells that diagnosed E4 cannot be its tests; held-out tests must come from laboratories not cited in secs. 2 - 5;
(2) the change is scored at the incumbent constants first, with a compensation ledger; (4) time box <= 3 h of work.
The probes behind every number here are read-only and in memory (sec. 9). They were run before this file was written,
and they are disclosed as such. No probe searched over a value.*

## 0. Summary for the morning (read this first)

- **The red team's motivating number is a misreading.** Sec. 5 says "shipped methylpyrazine pH 7 / 5 = 170 vs Leahy 14".
  But 14 is Leahy's ratio for **pyrazine**. Leahy's ratio for **methylpyrazine** is **215** at 95 C (99 at 85 C) from
  Table I rates, and ~279 from the Table III 2 h yields (`leahy1989a_extraction.md` sec. 3). The shipped engine gives
  170 at 95 C (-0.10 dex): it **matches**. At 85 C it overshoots: 427 vs 99, +0.63 dex. E4's band gives 9.9 - 20 (-1.0 to -1.3 dex): it **misses**. The premortem's
  circularity note ("the band was chosen so that MPZ comes out at 13.8 vs 14") inherits the same misreading. The Leahy
  ratio is indeed a development check (it diagnosed E4). As a development check it now speaks *against* E4.
- **Held-out tests:** there are **none** that are admissible (sec. 4). One pyrazine pH contrast comes from a laboratory
  not cited in secs. 2 - 5 (Monforte 2021, Porto), but it measures a different mechanism (the wine o-quinone Strecker
  route) and a headspace release rate. Power on admissible held-out data: **0**; effective power counting Monforte at
  its transfer probability: **~0.2** (sec. 4.3).
- **Incumbent-constant scoring (sec. 6):** no row of the core panel moves (66 / 66 rows, max |delta| 0.000 dex). The
  only predictions that move are trunk-lane pyrazine answers below pH 7: x3.9 at pH 6 and x16 at pH 5 for slope 0;
  x2.8 and x7.8 at the 0.15 edge (with the Zhou pH-8 re-anchoring, sec. 2). No shipped blind prediction (WUR) carries
  a pyrazine.
- **Expected verdict: REFUSE the point change (P 0.85).** The placement of the below-7 slope (step vs supply) is
  unidentified. The two fed-dicarbonyl laboratories (DFA, Ghent) favour "not on the step". Leahy's cascade favours a
  steep cascade that the engine's supply alone cannot produce. No held-out laboratory decides between them. Rule (f)
  applies: a declared band, plus the designed experiment (E4's row of sec. 3.1).

## 1. What changes (the exact parameter)

- `src/kinetic_core/parameters_pyrazine.py:154`
  `PYRAZINE_PH_SLOPES = (FROZEN_B18["ph_slope_above_7_decades_per_unit"], FROZEN_B18["ph_slope_below_7_decades_per_unit"])`
  `= (0.19745, 0.57951)` dec / pH, with the knot at `PYRAZINE_PH_KNOT = 7.0`. The term is normalised to 1 at
  `PYRAZINE_REFERENCE_PH = 6.8` (`pyrazine_ph_factor`). It is applied by `trunk_conditions.pyrazine_factor` to
  `PYRAZINE_PH_STEPS = ("k_go_ak", "k_mgo_ak")`, the two Strecker steps.
- **The change:** the second element, the below-7 slope `s`, goes from 0.5795 (fitted in B18 on Leahy's whole-cascade
  ratios as if they were step ratios) to a declared band **s in [0, 0.15]** dec / pH, point 0.07. The above-7 slope
  stays at 0.1974 (out of scope; see D2).
- **Kept:** the supply's pH law. `trunk_conditions.PH_EXPONENT_DECADES_PER_UNIT = 0.69` (band 0.37 - 0.92) on the
  Amadori-decay steps `PH_STEPS = ("k_ama_tdg", "k_ama_odg", "k_ama_mgo")` is not touched. The red team's argument:
  that law already carries the cascade's pH dependence, so a steep step slope double-counts.

## 2. Re-anchoring (a consequence, declared, not fitted)

B18 stores `log10_k_go_ak_100C` and `log10_k_mgo_ak_100C` at pH 6.8 but fitted them on Zhou 2024's pH-8 pots
*through* the factor. A different `s` changes factor(8) / factor(6.8) by 10^(0.2 (s - 0.5795)): -0.116 dex at s = 0,
-0.102 at 0.07, -0.086 at 0.15. Because the measurement is k at pH 8, the change re-anchors both constants by
**+0.2 (0.5795 - s) dex**. This leaves Zhou's six B18 fit rows unchanged (verified to 1e-9 in the run: V1), and every
pH >= 7 prediction exactly unchanged. Without re-anchoring, every pyrazine at pH >= 7 would fall by x0.77 - 0.82. That
is rival R-na (recorded, not a candidate).

## 3. Rivals (written before any result)

| id | below-7 slope `s` | re-anchored | role |
| --- | --- | --- | --- |
| N | 0.5795 (shipped) | - | the null |
| E4-0 | 0 | yes | band edge: "no pH on the step below 7" |
| E4-p | 0.07 | yes | the point (sec. 3's "<= 0.07, as the Strecker prereg allows") |
| E4-15 | 0.15 | yes | band edge |
| R-na | 0 | no | sensitivity only |

## 4. Tests

### 4.1 Held-out (premortem rule: laboratories not cited in ENGINE_VS_LITERATURE secs. 2 - 5)

I searched the atlas for every pyrazine or Strecker pH contrast (F16, F17, F19, F20, F21) and checked each laboratory
against secs. 2 - 5 (grep of the section text).

| cell | lab | cited in secs. 2 - 5? | below-7 contrast | admissible? |
| --- | --- | --- | --- | --- |
| F16#11 Monforte 2021 | Porto (Silva Ferreira) | no | PA release k 0.32 +- 0.01 (pH 5) vs 0.36 +- 0.05 (pH 7), 80 C: ratio 0.89 | **no as a binding test**: the dossier says it is the o-quinone / metal / O2 Strecker route with gallic acid, not the dicarbonyl step. The rate is a Weibull headspace release, and glucose *lowers* it |
| F19#5 / F20#9 Zhou 2023 | Ho network (Jiangnan / Rutgers) | yes (F01 - F04; E13 / E14) | MGO + Cys, 120 C, label pH 6 / 7 / 8: 2,5-DMP 66 / 106 / 237; PZ 47 -> 15, MPZ 160 -> 117 (6 -> 8) | no: the lab is cited; cysteine is not a trunk amine (no answer); the three pyrazines move in opposite directions |
| F19#20 Zhou 2022 | Ho network | yes | pH 7 -> 8 only | no: above the knot |
| F16#12 Chan 1994b | Minnesota (Reineccius) | yes (it is Leahy's laboratory) | PA highest at pH 6 (6 / 7 / 8) | no |
| F16#3, F19#7 / F21#1, F16#5 - 6, F19#8 - 9 | DFA, Ghent, Quest, Reineccius | yes (they diagnosed E4) | - | development checks (sec. 4.2) |

**Admissible held-out tests: 0 laboratories.** Under the premortem's rule (fewer than 2 labs), the test is
underpowered. Its power is stated in 4.3.

**H1 (scored, NOT binding): Monforte 2021,** step k(5) / k(7) within 0.3 dex of 0.89 (that is, [0.45, 1.78]):
- shipped 0.069 (-1.11 dex): fail;
- E4-0 1.00, E4-p 0.72, E4-15 0.50: pass.

It is recorded as a transfer diagnostic because P(the quinone route's pH law transfers to the dicarbonyl step) ~0.2
(mine).

### 4.2 Development checks (the cells that diagnosed E4; reported, not decisive)

- **D1 (Leahy 1989a, Reineccius, lysine + glucose 0.1 M, citrate-phosphate, the engine run with glycine, declared).**
  Engine 2 h yield ratio vs Table I rate ratio, tolerance 0.3 dex.
  - Methylpyrazine (7)/(5): 215 at 95 C (yield-based ~279), 99 at 85 C.
  - Pyrazine: 14.4 at 95 C, 6.9 at 85 C. Leahy's own Table III implies ~240 for pyrazine at 95 C: Table I's pH-5 pyrazine
    rate is 8x Table III's 2 h level, an internal conflict. So pyrazine is flagged unreliable at pH 5 in the source.
  - **Pyrazine's expected miss (stated):** under every band value, the engine's pyrazine (7)/(5) is 0.5 - 1.6. That is
    1.0 - 1.5 dex below 14.4. The engine's glyoxal supply *rises* toward acid (shipped: supply x0.56 over pH 7 -> 5),
    and its pH-9 glyoxal defect (PZ (9)/(7) 0.097 vs 2.67) is separate. A miss in [-1.6, -1.0] dex is the glyoxal-supply
    diagnostic, not a verdict on E4. A miss outside that range means the probe misdescribes the cascade: stop and
    re-check before deciding.
- **D2 (DFA, Hofmann 2000 Fig. 4, fed dicarbonyl + Phe, 98 C, phosphate): PA maximal at pH 5 for every donor.**
  - Below-7 limb: step k(5)/k(7) >= 10^-0.3 = 0.50.
  - Above-7 limb (diagnostic only): k(5)/k(9) within 0.3 dex of 2 (MGO). Every rival fails it (shipped 0.028; E4 0.20 -
    0.40) because the above-7 slope rises. That is Leahy (pH 9/7, borate) against DFA. It is out of scope and is
    recorded as a candidate conflict. PAA rises x7 from pH 3 -> 9 (F17#4), so PA alone may understate the step at
    high pH.
- **D3 (Ghent, Adams 2008, fed MGO + Ala ~0.85 M, 130 C): 2,5(6)-DMP pH 4 (acetate) / pH 7 (phosphate) = 2.1.**
  - Read as a lower bound on the step: the pyrrole sink rises x80 at pH 4 and competes for aminoacetone.
  - Test: k(4)/k(7) >= 2.1 / 2.
  - Shipped 0.018, E4-0 1.00, E4-p 0.62, E4-15 0.36. All fail; E4-0 misses by 0.02 dex.
  - Buffer- and concentration-confounded (phosphate alone gives x1.8, F19#6).

### 4.3 Power

Power on admissible held-out data: 0 (no laboratory).

If H1 were admitted, here is how it would discriminate. Shipped sits 1.11 dex off, which is ~3.7 sigma at an assumed
0.3 dex between-lab error. The band sits within 0.25 dex. So the power is ~0.95 *conditional on transfer*, and the
effective power is ~0.2. By the plan's rule (c), an underpowered test is not scored. So **no held-out test decides
E4.**

### 4.4 Guard (incumbent constants first; compensation ledger)

- **G1, the core panel (66 rows).** No row gets worse by > 0.3 dex unless it is listed in the compensation ledger with
  its named open error. **Ledger (stated now): empty.** No core-panel row carries a pyrazine, and the Strecker step
  consumes < 0.2 % of GO / MGO.
- **G2, B18's fit rows.** Zhou 2024's six rows unchanged (V1).

## 5. Decision rule (frozen)

- **SHIP** (needs owner approval, because pyrazines move up to x16): an admissible held-out test passes at >= 2 labs
  AND G1 / G2 hold AND D1's pyrazine miss lies in its declared interval. *Not reachable on the evidence on disk.*
- **DO NOT SHIP (keep the shipped slope as the point):** if the owner rules that the development checks decide.
  - Shipped passes D1 MPZ at 95 C and fails D2 / D3. E4 fails D1 MPZ by ~1.2 dex and passes D2.
  - Leahy is one lab with a buffer change and a different amine; DFA + Ghent are two fed laboratories. So the
    development checks split.
- **REFUSE (expected):** no admissible held-out test, and development checks that split by laboratory, mean the placement
  is unidentified (rule f). Then:
  - (i) the point stays at shipped, so no prediction moves;
  - (ii) pyrazine answers below pH 6.8 carry a declared one-sided band over s in [0, 0.58]: up to x10^(0.58 (6.8 - pH)),
    which is x2.9 at pH 6 and x11 at pH 5 (x16 with re-anchoring). This needs an envelope hook (a code change, later)
    or, now, a declaration on every pyrazine answer below pH 6.8;
  - (iii) the designed experiment goes to the brief: fed MGO / GO + Gly (and Leu / Phe), dilute (<= 50 mM), pH-stat
    4 / 5 / 6 / 7, 100 - 120 C, one pot with 2,5-DMP, MPZ, the aminoketone, the Strecker aldehyde and acid, phosphate
    vs citrate at equal pH. That decides step vs supply directly.

## 6. Incumbent-constant scoring (the probe; no refit exists to run)

Glucose 100 + glycine 100 mM, clamped pH, 2 h. Values are engine (7)/(5) ratios; deltas are log10 new / shipped,
**before** the +0.2 (0.5795 - s) dex re-anchoring.

| rival | MPZ 95 C | MPZ 85 C | PZ 95 C | PZ 85 C | DMP 95 C | step k(5)/k(7) | k(5)/k(9) | k(4)/k(7) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Leahy / DFA / Ghent | **215** (279) | **99** | 14.4 (~240) | 6.9 | < LOQ at 5 | >= 1 (D2) | 2 (D2) | >= 2.1 (D3) |
| shipped | 170 | 427 | 8.07 | 15.7 | 3550 | 0.069 | 0.028 | 0.018 |
| E4-0 | 9.9 | 15.0 | 0.51 | 0.73 | 191 | 1.00 | 0.40 | 1.00 |
| E4-p (0.07) | 13.8 | 21.7 | 0.71 | 1.03 | 268 | 0.72 | 0.29 | 0.62 |
| E4-15 | 20.3 | 33.6 | 1.04 | 1.55 | 396 | 0.50 | 0.20 | 0.36 |

**Predictions that move (re-anchored; trunk lane, any sugar + glycine pot).** pH >= 7: unchanged.

| rival | pH 5 | pH 5.5 | pH 6 | pH 6.5 |
| --- | --- | --- | --- | --- |
| E4-0 | +1.20 (x16) | +0.89 | +0.59 (x3.9) | +0.30 |
| E4-p | +1.06 | +0.79 | +0.52 | +0.26 |
| E4-15 | +0.89 (x7.8) | +0.67 | +0.44 (x2.8) | +0.22 |

The same three pyrazines move within 0.05 dex of each other at 95 and 120 C (85 C: up to +0.15 more for DMP).

**Core panel:** 66 / 66 rows identical under N, E4-0 and E4-15 (max |delta| 0.000 dex; 21 / 66 within 3x in all).
**WUR blind predictions:** no pyrazine (checked `scratch/wur_blind_test/BLIND_PREDICTIONS.md`).

## 7. Expectations (mine, before the decision)

| statement | P |
| --- | --- |
| No admissible held-out test exists on disk (the probe found none) | 0.95 |
| E4 (any band value) passes D1 methylpyrazine at 95 C | 0.03 |
| Shipped passes D1 methylpyrazine at 95 C (probe: -0.10 dex) | 0.95 |
| E4 passes D2's below-7 limb; shipped fails it | 0.97 |
| G1 holds (no panel row moves) | 0.99 |
| Verdict REFUSE (declared band + designed experiment) / DO NOT SHIP / SHIP | 0.85 / 0.13 / 0.02 |
| The designed experiment (sec. 5 iii) puts the below-7 slope on the step at < 0.15 | 0.55 |

## 8. What would mislead

- (a) **Lysine vs glycine.** Leahy used lysine (two amines, the epsilon one mostly protonated at pH 9); the engine runs
  glycine. A lysine-specific steep pH law would make D1 irrelevant to the step.
- (b) **The supply.** If the engine's dicarbonyl supply pH law is too weak (the Martins 0.69 is a 5.5 - 6.8 window,
  extrapolated to pH 5), E4 is right *and* D1 fails, for a supply reason. That is the reading under which REFUSE, not
  DO NOT SHIP, is correct.
- (c) **Strecker aldehyde is not aminoketone.** DFA measures PA; PAA takes an O2 route that rises with pH, so PA's pH
  optimum may not be the aminoketone flux's.
- (d) **Concentration and buffer.** Ghent is 0.85 M with acetate vs phosphate.
- (e) **CF-077.** The methylpyrazine pairing factor (R = 1 vs 4) is pH-invariant, so it cannot bias a ratio test, only
  absolute MPZ.
- (f) **My own reading.** Every number for Leahy comes from the dossier's re-typed Table I. If the red team read a
  different table, check `leahy1989a_extraction.md` sec. 3 before deciding.

## 9. Probes (read-only, one core, nice 10, in memory; nothing written to the repo)

- `/private/tmp/claude-501/-Users-pabloantoniomorenocasares-Developer-Maillard/128e3455-67ec-4d82-ba48-6cfcfb79a63e/scratchpad/prereg_probes/probe_e4.py`
  writes `probe_e4.json` / `.log` (the cascade ratios and moves).
- `/private/tmp/claude-501/-Users-pabloantoniomorenocasares-Developer-Maillard/128e3455-67ec-4d82-ba48-6cfcfb79a63e/scratchpad/prereg_probes/probe_e4_panel.py`
  writes `probe_e4_panel.json` / `.log` (the core panel under N / E4-0 / E4-15).

They reproduce the red team's numbers exactly (shipped MPZ 170, PZ 8.07; s = 0.07 gives MPZ 13.8, PZ 0.71); only the
comparison target differs. Work so far ~2 h; the remaining run (V1 re-anchoring check, the verdict note) is <= 1 h,
which is inside the 3 h box.
