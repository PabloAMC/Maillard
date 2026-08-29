# Kinetic core, Build Wave B2.4 — the DECLARED weighting, and an ensemble

Pre-registered in [`results/validation/kinetic_core_b2_4_prereg.md`](kinetic_core_b2_4_prereg.md), saved to disk before any fit ran. Every claim below is scored in §5 whether it held or not.

**This wave changes no stoichiometry, no species, no reaction, no benchmark, no engine input, no hold-out row and no pass band.** It makes an accidental weighting explicit, replaces a point estimate with an ensemble, adds a continuous statistic to each of the two scorecards, and repairs a broken citation.

## 1. The declared exchange rate

`E` = decades of level error that ONE pH unit of endpoint error is declared to be worth. Implemented as `sigma_ph = 0.35 / E` on the three `zhou_final_pH_*` rows and nowhere else.

| tag | E | σ_ph | basis |
|---|---:|---:|---|
| **W-SHIPPED** | 1.4 | 0.250 | B2.3's ACCIDENTAL rate made explicit: sigma_log 0.35 over sigma_ph 0.25 = 1.40 decades per pH unit. Reproduces B2.3's objective exactly and is the control -- it is what 'declaring the weighting' costs when the declared value is the one already in force. |
| **W-HALF** | 0.7 | 0.500 | Half the shipped rate, i.e. the three pH endpoints down-weighted by 4x IN COST. Chosen as the smallest intervention large enough to be visible against a 10-dof objective, and chosen for that reason and not for any property of its outcome. |
| **W-MEASURED** | 0.28 | 1.250 | The exchange rate the corpus measures ON THE FIT SIDE. Kumazawa 2003's six DECLARED FIT retention rungs move 2-furfurylthiol from 0.995 to 0.110 survival across pH 3.0-6.4: 0.956 decades over 3.4 pH units, i.e. 0.281 decades per pH unit. This prices a pH miss at the level consequence the fit corpus's own pH-response measurement says a pH miss has. It REPLACES a first draft of this value that was derived from Hofmann's hold-out pH ladder -- see the correction note in this module. |

> B2.3 never declared a rate. Its three pH rows sat at σ = 0.25 **pH units** in the same sum of squares as 55 rows at σ_log = 0.20–0.60 **log-folds**, and D1 §3 measured the consequence: one pH unit priced at ~9× a 3× level miss, with `zhou_final_pH_from_pH8` alone carrying 44% of the entire B2.2→B2.3 refit. W-SHIPPED reproduces that objective exactly and is the control.

## 2. The free set — 20 of 48, frozen 28

**20 free, 28 frozen at their B2.3 values.** The route is D1's own stated fallback and the reason is arithmetic: B2.3's multistart trace records 4731 s and 2280 s for two full-vector starts, and this container delivers about one core of real throughput.

| parameter | clause |
|---|---|
| `k_tdp_fur` | R1 identified |
| `k_nf_mft` | R1 identified |
| `k_nf_mp3p` | R1 identified |
| `k_mgo_mp` | R1 identified |
| `k_glc_ha` | R1 identified |
| `k_dimer_fft` | R1 identified |
| `k_fft_decay` | R1 identified |
| `k_fur_decay` | R1 identified |
| `k_pent_caramel` | R1 identified |
| `k_cys_thermal` | R1 identified |
| `k_thiolate_loss` | R1 identified |
| `ph_acid_yield_per_sink_event` | R2 pH-drift (price-taker of the weighting) |
| `ph_arp_secondary_ammonium_pKa` | R2 pH-drift (price-taker of the weighting) |
| `k_fur_fft` | R3 topology switch (moved >=3 decades B2.2->B2.3) |
| `k_osone_decay` | R3 topology switch (moved >=3 decades B2.2->B2.3) |
| `k_dimer_decay` | R3 topology switch (moved >=3 decades B2.2->B2.3) |
| `k_ttca_deg` | R3 topology switch (moved >=3 decades B2.2->B2.3) |
| `k_thiol_decay` | R3 topology switch (moved >=3 decades B2.2->B2.3) |
| `Ea_decay_thiol_sink` | R4 decay barrier (traded on/off its bound) |
| `Ea_decay_carbonyl_sink` | R4 decay barrier (traded on/off its bound) |

## 3. THE ENSEMBLE — the spread, not the best member

| weighting | members | converged | budget-exhausted | costs | **S = log₁₀(max/min)** | S, local arm | Σr²_level at best | G1 | G2 |
|---|---:|---:|---:|---|---:|---:|---:|---|---|
| **W-SHIPPED** | 6 | 2 | 4 | 8.15, 8.15, 11.11, 19.41, 15.60, 15.59 | **6.79e-08** | 6.79e-08 | 15.5 | yes | no |
| **W-HALF** | 6 | 5 | 1 | 7.82, 7.77, 7.77, 15.23, 15.23, 15.23 | **0.2924** | 0.002948 | 15.12 | yes | yes |
| **W-MEASURED** | 6 | 4 | 2 | 7.40, 7.40, 7.40, 14.90, 14.90, 14.98 | **0.304** | 0.000183 | 14.42 | yes | yes |

> **Total cost is not comparable across weightings** — they are three different objectives. `Σr²_level`, the sum of squared residuals over the 55 non-pH rows at their unchanged σ_log, is, and it is the only cross-weighting comparison this report makes.

## 4. THE SHIPPING CHOICE

**Shipped: W-HALF.**

> G1: all three zhou_final_pH_* within 1.0 pH unit at the ensemble-best member. G2: >=4 of 6 members converged. S = log10(max cost / min cost) over converged members; ship the qualifying weighting with the smallest S; ties within 0.05 break toward the largest E. Declared in results/validation/kinetic_core_b2_4_prereg.md sec. 3 before any fit ran. THE EXAM AND THE PANEL DO NOT ENTER IT.

Qualifying: W-HALF, W-MEASURED. Tie broken toward the largest E: True. Fallback branch used: False.

## 5. THE PRE-REGISTRATION, SCORED

**7 held, 3 falsified, 0 unscoreable.**

- **FALSIFIED** — W-1 spread at W-SHIPPED > 0.30 (log10 max/min cost, converged)
  - pooled S = 6.79e-08 over 2 converged of 6; local arm alone S = 6.79e-08; 4 budget-exhausted. Costs: [8.146, 8.146, 11.106, 19.414, 15.595, 15.595]. TWO GLOBAL STARTS 2005 AND 1942 COST-UNITS APART LANDED WITHIN 0.0003 OF EACH OTHER, so within the ~20 identified coordinates the objective really is single-basin; the B2.2->B2.3 scatter lives in the ~28 coordinates this wave FROZE and says nothing about. THE FALSIFICATION IS QUALIFIED, AND THE QUALIFICATION IS NOT IN THE CONTROL'S FAVOUR: at W-SHIPPED only 2 members survived A1's budget to enter S, so the statistic is a two-member statistic. The same seeds at the other two weightings DID converge and DID resolve a second basin -- W-HALF: 5 converged, S = 0.2924 (local arm alone 0.002948), costs [7.82, 7.77, 7.77, 15.23, 15.23, 15.23]; -- W-MEASURED: 4 converged, S = 0.304 (local arm alone 0.000183), costs [7.4, 7.4, 7.4, 14.9, 14.9, 14.98]. Read together: the near-basin structure is genuinely flat (local-arm S is 0.003 and 0.0002), but the GLOBAL arm finds a distinct attractor at about twice the cost at every weighting where it is allowed to converge. W-1 is falsified as written and the objective is still multi-modal over the free set.
- **HELD** — W-2 exam geometric-mean fold spans >= 1.5x across the W-SHIPPED ensemble
  - 6 members re-sat: [28.53, 36.39, 37.2, 37.22, 48.81, 48.84]; range 1.712x
- **FALSIFIED** — W-3 down-weighting moves k_fur_fft below -2.0 AND the FFT pH slope below 3.0 decades
  - W-HALF: k_fur_fft -0.5576 (needs < -2.0), pH3->pH7 FFT slope 4.892 decades (needs < 3.0); W-MEASURED: k_fur_fft -0.574 (needs < -2.0), pH3->pH7 FFT slope 4.906 decades (needs < 3.0)
- **HELD ON THE CORRECT REFERENCE, FALSIFIED AS WRITTEN** — W-4 the charge-conservation fix stays free (B2.2's parameters under B2.3's stoichiometry; claim as written: within 10% of 26.8x)
  - On D1's own 23-point pool the as-was geometric mean is 23.97x. THE CLAIM MIS-NAMED ITS REFERENCE: 26.8x is D1's `B22_ref` cell (B2.2 stoichiometry), while the cell the claim describes is D1's `S_only` (B2.3 stoichiometry), published at 23.97x. Against 23.97x this reproduces D1 to four significant figures; against the 26.8x actually written it is -10.5%, just outside the 10% band. Either way the substantive claim -- the conservation fix is free -- stands, and the as-was PAIRED MEDIAN is 10.63x against D1's 10.63x. On the CURRENT 27-point pool the same number is 33.94x; the pool grew because concurrent Wave B6 gave the core a lipid lane that answers four points D1's core declined.
- **HELD** — W-5 the Hofmann pH-3/pH-7 block never reaches 3/4 and the ladder stays humped -- re-weighting does NOT fix a shape defect
  - panel block: W-SHIPPED 0/4, W-HALF 0/4, W-MEASURED 0/4 -- never reaches 3/4. Three-rung ladder (pH 5 fit-side-anchored, pH 3 and pH 7 hold-out, scored by the exam): monotone falling at ['shipped/FFT', 'half/FFT', 'measured/FFT']; monotone on BOTH compounds -- which is what the pre-registration declares as the falsifier -- at no weighting. REPORTED AGAINST THE CLAIM'S OWN PROSE, NOT ONLY ITS FALSIFIER: the claim asserts the ladder 'stays humped ... at every weighting', and only MFT does. FFT is monotone FALLING at all three weightings, at 4.88-4.91 decades against a measured 1.28 -- the right SHAPE with the wrong STEEPNESS, roughly four decades too steep. So W-5's verdict holds and its substance holds -- re-weighting fixed neither the block nor the shape -- but its description of the defect is half wrong, and the FFT defect is a slope defect rather than a hump. LIMITATION: three rungs, not D1's seven-rung sweep, so a hump BETWEEN pH 3 and pH 5 is invisible to this test.
- **HELD** — W-6 fed_c2c3_MFT_pH3 is the largest single residual at every weighting
  - W-SHIPPED: fed_c2c3_MFT_pH3 at 2.087 sigma; W-HALF: fed_c2c3_MFT_pH3 at 2.087 sigma; W-MEASURED: fed_c2c3_MFT_pH3 at 2.087 sigma
- **HELD** — W-7 exam and panel median |log10 fold| agree in SIGN at every weighting
  - W-HALF: exam +0.000, panel -0.072 (same direction); W-MEASURED: exam +0.000, panel -0.110 (same direction). THE VERDICT IS DEGENERATE AND IS NOT EVIDENCE OF AGREEMENT: the exam's median |log10 fold| is BIT-IDENTICAL across all three weightings, so the exam contributes no sign and the test cannot fail on this statistic.. ON THE GEOMETRIC MEAN -- the other continuous statistic this wave added, and the one that actually moves -- THE TWO SCORECARDS DISAGREE IN SIGN at W-HALF, W-MEASURED: W-HALF exam +0.328x, panel -0.039x; W-MEASURED exam +0.765x, panel -0.562x. The panel improves as E falls while the exam worsens. By W-7's own stated consequence -- 'if they ever disagree in sign, the four shared rows are masking something and the wave says so' -- THE WAVE SAYS SO.
- **HELD** — W-8 sum_r2_level improves monotonically as E falls (W-SHIPPED >= W-HALF >= W-MEASURED)
  - W-SHIPPED: 15.5; W-HALF: 15.12; W-MEASURED: 14.42
- **FALSIFIED** — W-9 the pH-endpoint misses grow monotonically as E falls, AND W-MEASURED fails gate G1
  - W-SHIPPED: worst |miss| 0.1698 pH units; W-HALF: worst |miss| 0.2807 pH units; W-MEASURED: worst |miss| 0.6625 pH units. W-MEASURED gate G1 = PASS
- **HELD** — W-10 the pre-declared criterion does NOT select the best-exam weighting
  - criterion ships W-HALF; best exam geometric mean is W-SHIPPED at 48.84x. All: W-SHIPPED 48.84x, W-HALF 49.17x, W-MEASURED 49.61x
- **HELD** — PREREG sec.5 the rebuilt probe reproduces the DIRECTION of the B2.3 table (carry-both leaves all three Zhou endpoints above pH 4.5)
  - three encodings distinct: True; reproduces the B2.3 prereg table cell-by-cell: True; mean |miss| on the two acidifying anchors -- shipped 0.3407, carry-both 1.755 pH units

## 6. THE AMINE-FATE PROBE, REBUILT

D1 §7 found that `ph_state.AMINE_FATE_BASIS` cited a probe that **cannot run** — `KeyError: 'AMN'`, a species removed from the tree — and that two of its three encodings had collapsed onto one code path.

`scripts/generators/probe_amine_fate_b2_4.py` rebuilds the axis against the current species set, deriving the released amino nitrogen as `(Cys + ARP + TTCA at t=0) − (Cys + ARP + TTCA at t)` and adding it back as ammonium at pK_a 9.25. No new species, no new parameter, no reaction changed.

- three encodings genuinely distinct: **True**
- reproduces the B2.3 pre-registration's published table, cell by cell to two decimals: **True**
- mean |miss| on the two acidifying anchors (Zhou pH 6, pH 7): shipped **0.3407** pH units vs carry-both **1.755**

> **The defect was in the script, not in the evidence.** The B2.3 table's digits are reproducible from the current tree, so the declaration keeps its basis rather than being weakened; what changes is that the citation now points at a probe that runs, and the self-referential last sentence of `AMINE_FATE_BASIS` is gone.
