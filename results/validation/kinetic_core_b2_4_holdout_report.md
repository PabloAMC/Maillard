# Kinetic core, Build Wave B2.4 — the blind re-sit, under THREE declared weightings

Both scorecards, run once per declared weighting from the frozen ensemble-best parameters of that weighting. **No optimiser runs in either scorer, and neither scorer was forked**: the existing `generate_kinetic_core_b2_3_holdout` and `generate_cutover_final_exam` are called with two things rebound — which fit report the parameters come from, and which basename the artifact is written under.

## 1. THE FOUR SHARED ROWS — declared in both scorecards from this wave on

D1 §5 established that four of the exam's 23 answered points are **the same measurements** as four hold-out panel rows, not analogues of them. The panel reads `hofmann_ribose_pH7_FFT` at 0.0020× and the exam reads the same measurement at 499×: one number, inverted. **The two scorecards are therefore not independent evidence on this axis**, and agreement between them here is one measurement counted twice.

| panel row | exam point |
|---|---|
| `hofmann_ribose_pH3_FFT` | `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3` / 2-Furfurylthiol (FFT) |
| `hofmann_ribose_pH3_MFT` | `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH3` / 2-Methyl-3-furanthiol (MFT) |
| `hofmann_ribose_pH7_FFT` | `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7` / 2-Furfurylthiol (FFT) |
| `hofmann_ribose_pH7_MFT` | `mp_holdout_hofmann1998_ribose_cysteine_145C_20min_pH7` / 2-Methyl-3-furanthiol (MFT) |

- W-SHIPPED: panel declares 4 shared rows, exam declares 4.
- W-HALF: panel declares 4 shared rows, exam declares 4.
- W-MEASURED: panel declares 4 shared rows, exam declares 4.

## 2. THE WEIGHTING TABLE — fit, panel and exam under each declared value

| | E | fit cost | Σr²_level | **panel gating** | **panel median \|log₁₀\|** | panel geo-mean | **exam geo-mean (D1 pool, n=23)** | exam geo-mean (all 27) | exam paired median | exam in band | Hofmann pH block |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **W-SHIPPED** | 1.4 | 8.146 | 15.5 | **11/32** | **0.6909** | 13.38x | **48.84x** | 62.23x | 50.13x | 4/27 | 0/4 |
| **W-HALF** | 0.7 | 7.769 | 15.12 | **12/32** | **0.6185** | 13.34x | **49.17x** | 62.59x | 50.13x | 4/27 | 0/4 |
| **W-MEASURED** | 0.28 | 7.4 | 14.42 | **12/32** | **0.5807** | 12.81x | **49.61x** | 63.06x | 50.13x | 4/27 | 0/4 |

> **THE POOL MOVED, AND IT IS NOT THIS WAVE'S DOING.** D1's geometric means are over **23** answered points. The exam now answers **27**: concurrent Build Wave B6 gave the core a lipid lane, so four `matrix_path` points it used to decline are now scored — at a family median near 1900×. Every cross-wave comparison in this report is therefore made on **D1's own 23-point pool** (all scored points except the lipid family), with the 27-point number printed beside it. Comparing a B2.4 27-point mean against D1's 26.8× would be comparing two different pools and calling the difference a result.

> **Read the two continuous columns, not the two censored ones.** The panel's gating count is a pass/fail tally that cannot see a 100× degradation on a row that was already failing — D1 §5 measured 1.42 net decades of B2.2→B2.3 degradation that cost it nothing. The exam's paired median is a rank statistic over a lumpy 23-point pool that swung 4.7× between those same waves while the geometric mean moved 1.78×. Both continuous statistics are new in this wave and both are permanent.

## 3. THE HOFMANN pH LADDER — the shape defect, scored

Measured: FFT 229.0 / 121.0 / 12.0 ppb and MFT 553.0 / 198.0 / 25.0 ppb at pH 3 / 5 / 7 — **monotone falling on both**, 1.281 and 1.345 decades. The pH-5 rung is a declared FIT anchor; pH 3 and pH 7 are hold-out and appear here only because the exam has already scored them.

| weighting | cmpd | pred pH 3 | pred pH 5 | pred pH 7 | decades | monotone falling? |
|---|---|---:|---:|---:|---:|---|
| W-SHIPPED | FFT | 1581 | 102.3 | 0.02097 | 4.877 | yes |
| W-SHIPPED | MFT | 13.27 | 172.5 | 0.04587 | 2.461 | no |
| W-HALF | FFT | 1598 | 103.1 | 0.02051 | 4.892 | yes |
| W-HALF | MFT | 13.27 | 172.6 | 0.04626 | 2.458 | no |
| W-MEASURED | FFT | 1584 | 103.2 | 0.01967 | 4.906 | yes |
| W-MEASURED | MFT | 13.31 | 173.5 | 0.0463 | 2.459 | no |

> **LIMITATION, stated:** three rungs, not D1 §4's seven-rung sweep. A hump lying between pH 3 and pH 5 is invisible to this test. The three rungs used are the two the exam already scores and the one the fit already anchors, so nothing new was integrated at a hold-out condition to build this table.

## 4. THE ENSEMBLE ON THE EXAM — W-2's re-sits

| member | **exam geo-mean (D1 pool)** | exam geo-mean (all) | exam paired median | in band |
|---|---:|---:|---:|---:|
| `shipped_s0` | **48.84x** | 62.23x | 50.13x | 4/27 |
| `shipped_s1` | **48.81x** | 62.2x | 50.13x | 4/27 |
| `shipped_s2` | **36.39x** | 48.43x | 50.13x | 2/27 |
| `shipped_s3` | **28.53x** | 39.37x | 50.13x | 5/27 |
| `shipped_s4` | **37.2x** | 49.35x | 50.13x | 3/27 |
| `shipped_s5` | **37.22x** | 49.37x | 50.13x | 3/27 |

> This table is the wave's sharpest single result if the numbers spread: it says the exam score is a property of which basin the optimiser happened to find, not of the model.

## 5. THE PRE-REGISTRATION, SCORED

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
