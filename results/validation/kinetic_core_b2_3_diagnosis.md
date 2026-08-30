# Build Wave B2.3 — what the numbers said, AFTER the scorecards were frozen

**Everything in this file was written after looking at the fit, the hold-out
scorecard and the exam. NO PARAMETER WAS CHANGED, before or after.** The frozen
fit is `kinetic_core_b2_3_fit_report.json`, the scorecard is
`kinetic_core_b2_3_holdout_report.md`, the exam is `cutover_final_exam.md`, and
the pre-registration is `kinetic_core_b2_3_prereg.md` — written and saved to
disk before the refit ran. This file exists so that a post-hoc explanation is
filed as a post-hoc explanation rather than folded back into a report where it
would read like a prediction.

---

## 1. The headline

| | B2.1 | B2.2 | **B2.3** |
|---|---:|---:|---:|
| fit reduced chi-square | 2.01 | 1.80 | **1.635** |
| fit rows / free parameters | 51 / 44 | 58 / 48 | **58 / 48** |
| individually identified parameters | — | 3 / 48 | **11 / 48** |
| gating hold-out (B2.2 denominator) | **15 / 33** | **12 / 33** | **12 / 33** |
| gating hold-out after the permanent demotion | — | — | **11 / 32** |
| rows FIXED / REGRESSED vs predecessor | 11 / 0 | 6 / 9 | **3 / 4** |
| Kang cooled endpoint, 100 / 120 °C (measured 4.9) | — | **11.42 / 11.42** | **5.60 / 5.33** |
| Zhou drift anchors, worst error (3 declared FIT) | — | 0.08 | **0.16** |
| `arp_secondary_ammonium_pKa` | — | 10.77 | **7.50** |
| exam, core paired median fold | 24.93× | **10.65×** | **50.13×** |
| exam, old-lane paired median | 12.65× | 12.65× | 12.65× |
| exam, Yiltirak family median (completed / as-was) | 290.5× | 262.3× | **193.5× / 444.5×** |
| exam, worst single point | 1475× | 1475× | **648× / 2380×** |

**The wave did exactly the thing it was sent to do and paid for it on the
exam.** The charge leak is gone — Kang's pot moved from six pH units wrong to
about half a unit — the Amadori pKa fell off its unphysical perch onto the value
the chemistry demands, three times as many parameters became identified, and the
fit improved. Against that, the exam's paired median went from 10.65× to 50.13×
and the core no longer beats the old lane on the only apples-to-apples number in
that report. **Neither half of that should be reported without the other.**

---

## 2. Pre-registration scorecard — every claim, held or not

**8 held, 5 falsified, 4 half-falsified, 3 not-yet-decidable.**

### The fit

| # | claim | outcome |
|---|---|---|
| F-1 | χ²ᵣ between 1.5 and 4.0, and **worse** than B2.2's 1.80 | **HALF-FALSIFIED.** 1.635 is inside the band but it is *better*, not worse. The same half-falsification B2.2 recorded on its own F-1, and for the same reason: a wave that removes a defect keeps expecting to pay for it and keeps not paying. |
| F-2 | the three Zhou drift endpoints reproduced within 0.5 pH units | **HELD.** Worst error **0.16 units** (pH-7 anchor: 3.584 predicted against 3.42). Predicted 3.076 / 3.584 / 5.054 against 3.22 / 3.42 / 5.07. Worse than B2.2's 0.08 — the fix costs a little precision on its own calibration target, which is what carrying a real constraint instead of a free nuisance parameter looks like. |
| F-3 | `arp_secondary_ammonium_pKa` falls **below 9.5**, toward the 7–8 an N-glycosylated secondary amine should give | **HELD, and this is the wave's most informative single number.** **7.50**, from B2.2's 10.77. See §3. |
| F-4 | `acid_yield` falls below 0.90 | **FALSIFIED.** 0.9675, essentially unmoved from B2.2's 0.9679. |
| F-5 | `Ea_decay_thiol_sink` stays above 200 kJ/mol, still bound-limited | **HELD**, 216.1 — but it is **no longer pressed against its 250 ceiling** (B2.2: 248.0, clearing the bound by 2 kJ/mol). It has become an estimate rather than a bound. |
| F-6 | 2–8 of 48 parameters individually identified | **FALSIFIED, in the good direction. 11 of 48**, against B2.2's 3. See §4. |
| F-7 | both Kang pots' cooled endpoints below pH 7.0 | **HELD, with room to spare: 5.60 and 5.33**, against B2.2's 11.42, against a measured 4.9. |
| F-8 | the Hofmann 0.5 M phosphate pot still finishes below pH 4.5 cooled | **HELD.** 3.95. Buffers in this corpus still do not hold, and B2.2's §5 contradiction 1 stands unchanged. |

### Blind re-sit (a) — the full hold-out panel

| # | claim | outcome |
|---|---|---|
| H-1 | gating between 10 and 16 of 33 | **HELD.** 12 / 33 on B2.2's denominator — **exactly level with B2.2**. 11 / 32 after Amendment 9's permanent demotion. |
| H-2 | 3–12 rows change status | **HELD.** Seven: 3 FIXED, 4 REGRESSED. |
| H-3 | the Zhou pH-8 block improves in aggregate | **HELD.** Median \|log₁₀ fold\| **0.636 → 0.474**; mean **1.189 → 0.881**. `zhou_pH8_MFTD` 4.33× → 0.47× and `zhou_pH8_FUR_arp_alone` 7.46× → 0.19×. |
| H-4 | the `kang_switch_on_*` rows still fail | **HELD.** Both. |
| H-5 | `cerny2007_branch_responsiveness` still fails | **HELD** (0.300×). |
| H-6 | at most 2 of the 4 Zhou dimer rows pass | **HELD.** One (`zhou_120C_dimer_share_pH8`, 12.6× → 0.72×). The invariance row improved 72.4× → 11.0× and still fails. |
| H-7 | `vanseeventer_50C_MFT_per_day` and `hofmann2002_brew_80C_FFT` both still pass | **HELD.** 1.68× and 2.79×. |

### Blind re-sit (b) — the cutover exam, both ways

| # | claim | outcome |
|---|---|---|
| X-1 | the buffer-completed Yiltirak median is ≥ 3× better than as-was | **FALSIFIED, narrowly. 2.30×** (444.5× → 193.5×). Every one of the eight rows improved; the size was over-claimed, not the direction. |
| X-2 | Yiltirak stays 0 / 8 buffer-completed | **HELD.** |
| X-3 | the as-was column is not uniformly worse | **HELD, and more strongly than expected.** Four of the seventeen moved rows are WORSE buffer-completed, all Hofmann, and the overall as-was median (47.1×) is *better* than the completed one (50.1×). See §5. |
| X-4 | the core's paired median stays in 5×–25× in both columns | **FALSIFIED.** 50.1× completed, 47.1× as-was. |
| X-5 | the core still beats the old lane on the paired median | **FALSIFIED.** 50.1× against the old lane's 12.65×. B2.2's win is given back. |
| X-6 | the 17 declensions stay 17 | **HELD.** |
| X-7 | the Chang two-arm defect is not fixed | **HELD.** Both arms still predict 4041.10 ppb at 30 min, against measured 1459 (acetate) and 832 (water), even though the two bundles now carry *different* buffer blocks. The acrylamide lane has no pH term, so the completed field cannot reach it. |

### Blind re-sit (c) — Kang 140 °C, scored and NOT counted

| # | claim | outcome |
|---|---|---|
| K-1 | the 140 °C pH trajectory is no longer six units wrong | **HELD.** Cooled endpoint **5.24**, against B2.2's 5.42 at this rung and 11.42 at the two rungs that were actually fitted. |
| K-2 | `kang_140C_MFT` **degrades** to 2×–30× | **FALSIFIED, and it was the second-sharpest prediction in the file.** It **improved** to 0.467× from 2.15×, and would have passed again. See §6. |
| K-3 | `kang_140C_FFT` stays failing | **HELD** (0.0487×). |
| K-4 | `kang_140C_cys_conversion` stays passing | **HELD** (1.069×). |

### Buffer-field completion

| # | claim | outcome |
|---|---|---|
| B-1 | at least 6 of 21 bundles resolve to `buffer_unknown` | **FALSIFIED. Two** (Li 2026, Steinhagen 2021). The estimate was made from "only 2 of 8 source papers are on disk" and ignored that the repo already carries verbatim methods quotations for most of the rest. |
| B-2 | at least one paper affirmatively states NO buffer | **HELD. Seven** do. |
| B-3 | every bundle byte-identical outside `conditions.buffer` | **HELD**, proven by `tests/unit/test_kinetic_core_b2_3.py::test_the_buffer_completion_changed_nothing_but_the_buffer` against hashes taken before the completion script first ran. |
| B-4 | at least one bundle's `conditions.ph` disagrees with its source | **HELD. Five.** See §7. |

---

## 3. THE AMADORI pKa FELL FROM 10.77 TO 7.50, AND THAT IS THE WAVE'S RESULT

B2.2's diagnosis §5 contradiction 2 said this, before B2.3 existed:

> The calibrated Amadori ammonium pKa (10.77) is not chemistry. An
> N-glycosylated secondary amine should sit near 7–8. The fit has an incentive
> to push it high because that makes the Amadori electrically invisible and lets
> the strong-ion difference be set by cysteine alone. **It is a fitted nuisance
> parameter that happens to live in a pKa slot, and it must never be quoted as
> an acid-base constant.**

That was a *diagnosis of a mechanism*, and it made a prediction it could not
test: remove the leak and the incentive disappears. **B2.3 removed the leak and
the pKa went to 7.50** — inside the 7–8 band the contradiction named, without
anything in this wave pushing it there and without its bounds (5.0–11.0)
changing.

The mechanism is now legible. B2.2's four Amadori enolisations deleted the
Amadori's own amino-acid moiety, carboxylate and all, into untitratable
fragment pools. On Zhou's pot — whose entire charge inventory is 20 mmol/L of
fed ARP — that is the *dominant* leak, and it is a leak **diagnosis §3 did not
name**: §3 listed the four cysteine sinks, which at pH 5–8 are nearly
charge-neutral because an α-amino acid is a zwitterion. The optimiser could not
repair a stoichiometry, so it did the only thing available to it: it raised the
pKa until the Amadori carried no charge to lose.

**This is now a defensible acid-base constant rather than a nuisance value**, and
the qualification B2.2 attached to it can be lifted — with the caveat that it is
still *fitted on three endpoints*, not measured.

---

## 4. ELEVEN PARAMETERS IDENTIFIED, AGAINST THREE

`k_tdp_fur` · `k_nf_mft` · `k_nf_mp3p` · `k_mgo_mp` · `k_glc_ha` ·
`k_dimer_fft` · `k_fft_decay` · `k_fur_decay` · `k_pent_caramel` ·
`k_cys_thermal` · `k_thiolate_loss`

On 58 rows with 10 degrees of freedom, 11 of 48 is still a badly
under-determined fit and none of these should be lifted into another module as a
measurement. But the direction is informative and it was not predicted: removing
a bookkeeping degree of freedom that the optimiser had been *using* — the pKa —
sharpens everything downstream of the pH. Both decay barriers and both drift
constants remain unidentified.

`Ea_decay_carbonyl_sink` **moved onto its ceiling** (150.2 → 249.92 against a
250 bound) as `Ea_decay_thiol_sink` came off it (248.0 → 216.1). One
bound-limited decay barrier has been traded for the other, which is not
progress; it is the same non-identification wearing a different hat, and it is
recorded here rather than reported as a shift in the chemistry.

---

## 5. THE BUFFER FIELD HELPED YILTIRAK AND HURT HOFMANN, AND BOTH ARE REAL

Seventeen of the 40 exam points moved; 13 got closer and 4 got further.

* **All eight Yiltirak rows improved**, by 1.23×–3.90×. The FFT rows carry it:
  2337× → 599×, 2380× → 648×, 1596× → 469×, 705× → 226×. The exam's worst single
  point falls from 2380× to 648×. **B2.2's §4 was right that a large part of the
  Yiltirak blow-up was a data-schema gap and not chemistry** — and X-2 is right
  that what remains is still chemistry: 0 / 8, median 193×.
* **Four Hofmann rows got worse**, three of them badly (31× → 318×, 35× → 499×,
  21× → 517×). Five others improved. The Hofmann family's *median* improves
  (17.2× → 12.2×) and its in-band count goes 1 → 2, so the family is better
  while three of its points are much worse.

The reason the two families respond in opposite directions is not mysterious and
it is worth stating, because it is the argument for keeping both columns
permanently. 0.5 M phosphate at pH 5 has a van Slyke capacity of only ≈8.8 mmol
L⁻¹ per pH unit against a modelled acid load in the tens; it does not hold the
pH, it *shifts the trajectory*. Whether that shift helps a given row depends on
which side of the measurement the as-was prediction sat. **Four of the moved
rows improved by cancellation of two errors before the buffer was declared, and
would have looked like successes.** An artifact that reported only the completed
column would have hidden that, which is exactly why Amendment 9 made the
two-column format permanent.

---

## 6. THE KANG 140 °C PREDICTION WAS WRONG, AND THE ROW STILL SHOULD NOT COUNT

K-2 pre-registered that repairing the trajectory would *cost* `kang_140C_MFT`,
on B2.2's own reasoning that its pass was "a coincidence of a broken
trajectory". The row went from 2.15× to **0.467×** — it improved, and it would
have passed again.

**The demotion stands regardless, and this is the point of making it permanent
in advance.** Amendment 9 demoted the row because the defect was found *through*
its scoring, not because the pass was expected to evaporate. A rule that only
binds when it costs nothing is not a rule. The row is scored, printed, and not
counted.

What the outcome does tell us is that B2.2's §3 explanation was **half right**.
Its claim that the 140 °C pot's pH was six units wrong was correct and is now
fixed. Its inference that the pass therefore depended on the broken trajectory
was wrong: the pass survives the repair, so the mechanism carrying it is not the
pH error. That remains unexplained and is not explained here.

---

## 7. FIVE BUNDLES' RECORDED pH DISAGREES WITH THEIR OWN SOURCE — reported, not acted on

Amendment 9 licenses completing the **condition record**. Editing an executable
condition is a different act: it changes a prediction. All five are written into
the bundles' new `conditions.buffer.ph_disagreement` field and left alone.

1. **`external_validation_bi_2020_raw_pea_hexanal`** and **`..._roasted_...`** —
   `ph: 6.0`, and **Bi 2020 states no pH anywhere** (every page of the PDF
   searched; zero hits). Both are unsourced assumptions.
2. **`external_validation_liu_2023_ppi_offnote_baseline`** — `ph: 6.0` against a
   source range of 6.3–7.3 (mean ≈6.8). Already self-declared in the protocol
   YAML and still open.
3. **`external_validation_li_2026_spi_wg_hme_control`** — `ph: 7.0` is the
   *wheat-gluten enzymatic pretreatment* pH at 30 °C, a different unit
   operation, transposed onto a 160 °C extrusion run.
4. **`mp_holdout_glucose_only_autoclave_121C_Steinhagen2021`** — **and this one
   breaks the scale convention.** `ph: 4.36` is the pH measured *after*
   autoclaving; the initial value is 4.98. Every other bundle's pH is an initial
   bench reading, so this bundle is on a different scale from the other twenty,
   and a model handed 4.36 as an initial condition is being given part of the
   answer to its own trajectory. Recorded with `initial_ph_scale:
   at_temperature` so the anomaly is machine-visible.

---

## 8. Contradictions found — reported, not adjudicated

1. **THE EXAM GOT MUCH WORSE WHILE THE PANEL STAYED LEVEL.** Gating hold-out is
   flat at 12/33; the exam's paired median went 10.65× → 50.13×. The two
   scorecards disagree about this wave, and nothing here reconciles them. The
   honest reading is that the exam is dominated by the Yiltirak and Hofmann
   families, where the buffer field and the pH trajectory both moved hardest,
   while the hold-out panel is spread over systems where neither did.
2. **The core's win over the old lane is given back.** B2.2's X-5 recorded the
   first time the core beat the old lane on the paired median. B2.3 loses it,
   50.13× against 12.65×. That is a real regression on the only apples-to-apples
   number in the exam and it is the first thing a reader of the exam should be
   told; `generate_cutover_final_exam.py`'s headline finding is now COMPUTED
   from the rows rather than typed, so it says so on its own.
3. **`Ea_decay_carbonyl_sink` is now the bound-limited one** (249.92 against a
   250 ceiling), having traded places with `Ea_decay_thiol_sink`. Neither is a
   barrier; both are bound artefacts of an under-determined objective.
4. **The liberated-amino-nitrogen declaration is the wave's largest untested
   assumption.** `ph_state.AMINE_FATE_BASIS` states it, evidences it on
   Amendment 7's own anchors, and names the direction of the error if it is
   wrong (pots modelled too acidic, worst where amino-nitrogen loading is
   highest). It is not a measurement and it must not be quoted as one.
5. **`H2S` is still not in the charge balance** (pKa₁ 7.05, present at mmol/L).
   `ph_state.UNTRACKED_TITRATABLE` records it as the largest remaining gap and
   declines to close it, because closing it is a modelling addition rather than
   a conservation fix and would need its own declaration. **Thiamine is the
   second** — its thiazolium is a permanent cation and Zhang's pot charges 44.5
   mmol/L of it.
6. **The trunk and acrylamide lanes now carry the audit but not the balance.**
   Glycine, the Schiff base, the Amadori compound and melanoidin nitrogen are all
   titratable and all invisible to the charge balance, because those lanes have
   no pH state. The audit passes there trivially. The day either lane gets a pH
   term, six declared gaps become six defects.

---

## 9. Process notes

* **Two starts, seeded, ftol 1e-6 — the settings Amendment 9 allows and the ones
  B2.2 actually shipped.** Start 0 converged to cost **8.175** in 4731 s and is
  the frozen point; start 1 to 10.381 in 2280 s. Both are reproducible from seed
  20260828. B2.2's two starts gave 16.59 and 8.99, so B2.3's *worse* start beats
  B2.2's worse start and its better start beats B2.2's better one.
* **The memory ceiling is still the binding constraint.** The container peaked at
  6.5 GiB of 7.8 GiB during start 1 and the run took 117 min. A third start was
  not budgeted and could not have been run.
* **A pre-freeze probe was run before the amine-fate choice was made**
  (`scratch/b23_encoding_probe.py`), comparing three admissible encodings at
  B2.2's own frozen parameters against Amendment 7's three declared FIT anchors
  and against Kang's 4.9 out-of-sample diagnostic. It is disclosed in the
  pre-registration §0 rather than here, because it happened before the choice
  and not after it.
* **Two scorer defects B2.2 reported and could not fix under its
  pure-re-scoring mandate are fixed here** and the fixes both move labels in the
  conservative direction: the compound prereg claim now scores *both* halves
  (band and not-better-than-old-lane) instead of only the band, and the exam's
  narrative findings now compute their headline numbers from the rows instead of
  carrying B5's hard-coded 24.93× / 12.65× / 290×.
