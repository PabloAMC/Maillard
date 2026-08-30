# Diagnosis Wave D1 — reconciling the B2.3 exam regression with a level hold-out panel

**READ-ONLY WAVE. No file under `src/` or `data/` was modified, no parameter
moved, no benchmark was edited, and no refit was run.** Every number below comes
from `engine.predict(..., parameters=...)` with a `ProcessSpec.ph_drift`
override — the public API — or from evaluating the frozen B2.3 fit objective at
a point. Hold-out values were read only where the existing exam and hold-out
scorers already print them.

Machine-readable companion: `d1_exam_panel_reconciliation.json`.
Probes: `scratch/d1_factor_isolation.py`, `scratch/d1_objective_flatness.py`,
`scratch/d1_param_swap.py`, `scratch/d1_ph_response.py` (all untracked).

> **PROVENANCE, and a warning about the tree.** All four probes ran between
> 09:40 and 09:49 on 2026-08-29 against `audit-remediation` @ `7f65cca` with
> `src/kinetic_core/` untouched since 07:16 — `engine.py`, `ph_state.py`,
> `sulfur.py`, `parameters_sulfur.py` all at their B2.3 commit state. The
> reproduction check in §1 confirms it: the cells reproduce the published B2.2
> exam, the B2.3 as-was column and the B2.3 buffer-completed column exactly.
> **A concurrent Build Wave B6 began editing `src/kinetic_core/engine.py` and
> `__init__.py` at 09:53, adding a lipid-oxidation lane, i.e. after every number
> in this file was computed.** Nothing here is invalidated by that, but anyone
> re-running these probes on a later tree should expect a different import graph
> and should re-check the §1 reproduction row before trusting a comparison.

---

## 0. The answer in five lines

1. **The charge-declaration fix is exonerated.** At B2.2's own frozen
   parameters, swapping in B2.3's stoichiometry moves the exam's paired median
   from **10.65× to 10.63×**, leaves the in-band count at 5/23, and *improves*
   both family medians. Its share of the regression is **−0.1%**.
2. **The refit is 96% of it.** Holding stoichiometry and inputs fixed and
   swapping B2.2's parameter vector for B2.3's takes the paired median
   **10.63× → 47.12×**.
3. **The buffer-field completion is 4%** (47.12× → 50.13×) and *improves* both
   family medians while doing it — it is a rank effect, not an accuracy effect.
4. **The two scorecards never disagreed.** Four of the exam's 23 answered points
   are the *same measurements* as four hold-out panel rows, and both scorecards
   record the same degradation. The panel's gating count cannot see it: 22 of its
   34 gating rows were already failing, and **1.42 net decades of degradation
   landed entirely on already-failing rows**, where a pass/fail count is blind.
5. **The defect is in the objective, not in the chemistry.** The refit was
   *forced* — B2.2's vector scores χ²ᵣ **11.04** under the corrected
   stoichiometry — and 95% of the cost it recovered came from **six rows**, three
   of them pH endpoints weighted at σ = 0.25 pH units. It paid in a currency the
   objective does not price: the **formation-vs-pH shape**, which the fit has *no
   usable row for* and the hold-out and exam are both dominated by.

---

## 1. Factor isolation — the three-cell decomposition

All cells run the exam's own `_core_spec` builder over the nine sulfur bundles;
the five acrylamide points are carried unchanged so the paired median is over the
same 23 points the exam reports. **The reproduction check passes exactly:**
`B22_ref` reproduces every B2.2 exam number, `S+P` reproduces the B2.3 as-was
column, `S+P+B` reproduces the shipped column.

| cell | stoich | params | buffer | paired median | geo-mean fold | in band /23 | Hofmann med | Yiltirak med |
|---|---|---|---|---:|---:|---:|---:|---:|
| `B22_ref` | B2.2 | B2.2 | as-was | **10.65×** | 26.80× | 5 | 4.16× | 262.3× |
| `S_only` | **B2.3** | B2.2 | as-was | **10.63×** | 23.97× | 5 | 3.33× | 226.6× |
| `S+P` | B2.3 | **B2.3** | as-was | **47.12×** | 55.32× | 2 | 17.15× | 444.5× |
| `S+P+B` (shipped) | B2.3 | B2.3 | **completed** | **50.13×** | 47.67× | 3 | 12.21× | 193.5× |
| `B22_buf` | B2.2 | B2.2 | completed | 21.36× | 29.63× | 4 | 4.97× | 382.6× |
| `CBX_off` | carboxyl NOT carried | B2.3 | completed | 50.13× | 49.85× | 4 | 12.02× | 216.0× |

**Attribution, in decades of the paired median (total +0.673):**

| factor | decades | share |
|---|---:|---:|
| charge stoichiometry (`CBX` carried, centres declared) | **−0.001** | **−0.1%** |
| the refit (B2.2 → B2.3 parameter vector) | **+0.647** | **+96.1%** |
| buffer-field completion | **+0.027** | **+4.0%** |

Two things follow immediately, and they contradict the diagnosis's own §1
framing ("the wave did exactly the thing it was sent to do and paid for it on
the exam"):

* **The wave did not pay for the conservation fix.** `S_only` is the
  conservation fix in isolation and it is free — arguably slightly positive
  (geo-mean 26.8× → 24.0×, Yiltirak 262× → 227×, Hofmann 4.16× → 3.33×).
* **B2.2's 10.65× was NOT a compensation artifact of the charge leak.** That was
  the mandate's stated worry and the evidence says no: the number survives the
  conservation fix at 10.63×. What it does not survive is the refit.

The `CBX_off` cell says the same thing from the other side: at B2.3's
parameters, turning the carried carboxyl *off* changes the paired median not at
all (50.13× either way). The stoichiometry change that this wave is named after
has essentially no exam leverage in either direction.

---

## 2. Row level — the regression is concentrated in four rows, and it is a cliff

The Hofmann bundle names truncate to 46 characters in the exam's own table, which
hides that `mp_holdout_hofmann1998_{glucose,ribose}_cysteine_145C_20min_*` are
**two bundles each — pH 3 and pH 7**. That is the whole story.

| bundle | cmpd | declared pH | measured ppb | B2.2 | S_only | B2.3 as-was | B2.3 shipped | in-situ pH end |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| ribose | MFT | **7** | 25 | 22.4 (1.12×) | 53.4 | 1.18 | **0.0483 (517×)** | 6.70 |
| ribose | FFT | **7** | 12 | 1.59 (7.55×) | 4.79 | 0.342 | **0.0241 (499×)** | 6.70 |
| glucose | FFT | **7** | 6 | 28.0 (4.67×) | 27.5 | 187 | **0.0189 (318×)** | 6.59 |
| glucose | MFT | **7** | 4 | 35.9 (8.97×) | 35.1 | 35.3 | 10.6 (2.64×) | 6.59 |
| ribose | FFT | 3 | 229 | 77.5 (2.96×) | 75.8 | 1594 | 1660 (7.25×) | 3.37 |
| ribose | MFT | 3 | 553 | 5.99 (92×) | 5.99 | 2.67 | 13.1 (42×) | 3.37 |
| glucose | FFT | 3 | 7 | 25.6 (3.65×) | 25.4 | 420 | 94.2 (13.5×) | 3.31 |
| glucose | MFT | 3 | 3 | 31.9 (10.6×) | 31.9 | 32.3 | 32.9 (11.0×) | 3.31 |
| xylose | FFT | 5 | 96 | 91.2 (1.05×) | 89.2 | 1264 | 295 (3.08×) | 4.19 |
| xylose | MFT | 5 | 143 | 361 (2.53×) | 327 | 419 | 340 (2.38×) | 4.19 |
| Yiltirak ×8 | — | 5.5 | — | 262× med | 227× | 445× | 193× | 4.54–4.64 |

**Exactly four rows crossed the middle of the pool**, and the count of sub-10×
rows fell 11 → 7 with nothing entering. Three of the four crossings are the same
event: a prediction **collapsing to 0.02–0.05 ppb** — three to four orders of
magnitude — in every pot whose buffer-completed in-situ pH ends **above 6.5**.
This is a cliff, not a drift.

**The reported median row did not move, in either wave.** B2.2's 10.65× median
was `glucose pH3 / MFT`, which is 10.97× in B2.3. B2.3's 50.13× median is
`mp_holdout_glucose_asparagine_180C_Ye2024 / Acrylamide` — an **acrylamide-lane
row that this wave provably did not touch** (identical in all three columns; the
lane has no pH term). The paired median moved 4.7× because the rank order
re-sorted around it, not because a typical row got 4.7× worse. The geometric-mean
fold over the same 23 points moves **26.8× → 47.7×, a factor of 1.78** — still a
real regression, and the honest size of it.

---

## 3. Why the refit went where it went — it was forced, and by six rows

Evaluating (never refitting) **B2.3's own objective** at the two frozen vectors:

| | cost | χ²ᵣ (58 rows, 48 params, 10 dof) |
|---|---:|---:|
| B2.2's frozen x, under B2.3's stoichiometry | 55.21 | **11.04** |
| B2.3's frozen x | 8.18 | **1.635** |

B2.2's vector is *not* a near-optimum of B2.3's objective, so this was not a
lateral move across a flat landscape. But 95% of the Σr² it recovered (94.1
total) comes from six rows:

| row | r at B2.2's x | r at B2.3's x | Δr² | what it is |
|---|---:|---:|---:|---|
| `zhou_final_pH_from_pH8` | −6.430 | −0.065 | **41.34** | pH endpoint, **σ = 0.25 pH units** |
| `kang_120C_MFT` | 4.678 | 0.316 | 21.78 | one pot, one temperature |
| `kang_120C_FUR` | 3.308 | −0.345 | 10.82 | same pot |
| `kang_120C_FFT` | 2.902 | 0.222 | 8.37 | same pot |
| `zhou_pH7_dimer_over_MFT` | −2.545 | −1.483 | 4.28 | ratio |
| `zhou_final_pH_from_pH7` | −1.689 | 0.658 | 2.42 | pH endpoint, σ = 0.25 |
| | | | **89.0 of 94.1** | |

**One pH-endpoint row carries 44% of the entire refit.** The three
`zhou_final_pH_*` rows are scored in pH units at σ = 0.25, so a 1.6-unit miss is
6.4σ and costs r² = 41. A **10× fold error** on a Hofmann level anchor at
σ_log = 0.35 costs r² = 8.2. **One pH unit is priced at roughly nine times a 3×
level miss.** That weighting was never declared as a weighting; it is a
consequence of two rows being measured in different units with independently
chosen σ.

### What the refit paid with

The eight rows that got *worse* are dominated by the pH-shape rows and by the
Hofmann family the exam scores:

| row | r at B2.2's x | r at B2.3's x |
|---|---:|---:|
| `fed_c2c3_MFT_pH7` | −0.709 | **−1.241** |
| `fed_c2c3_MFT_pH3` | 1.875 | **2.087** ← largest residual in B2.3's whole fit |
| `vanseeventer_cys_conversion` | 0.432 | 0.678 |
| `whitfield_nf_cys_MFT` | −0.315 | −0.554 |
| `hofmann_glucose_FFT` | −0.177 | −0.331 |
| `hofmann_ribose_FFT` | 0.051 | −0.144 |
| `hofmann_fructose_FFT` | 0.228 | 0.330 |
| `hofmann_ribose_MFT` | 0.157 | −0.189 |

**Every Hofmann FFT anchor in the fit got worse, and both of the fit's only two
direct formation-vs-pH rows got worse.** And `generate_kinetic_core_b2_3_fit.py`
declares those two rows unfittable in its own note:

> `fed_c2c3_MFT_pH7` … **"THE C2+C3 ROUTE HAS NO pH FACTOR IN THIS MODULE.
> Assigning one would be inventing a mechanism to fit two numbers, so these two
> rows are carried in the objective and the residual is REPORTED."**

So the objective contains **zero usable constraint on how thiol *formation*
responds to pH.** What it does constrain on the pH axis is (a) thiol
*destruction*, through Kumazawa's six retention rungs (pH 3.0–6.4, 121 °C,
10 min), and (b) the pH *trajectory*, through Zhou's three endpoints. Every
Hofmann anchor in the fit is **pH 5.0 only**. The exam's Hofmann hold-out rows
are pH **3 and 7** — the two ends of an axis the fit never looked at.

### How far the parameters moved for a 10% cost gain

| constant | B2.2 | B2.3 | move |
|---|---:|---:|---|
| `k_osone_decay` | −0.858 | −8.696 | **÷ 7×10⁷** (channel switched off) |
| `k_dimer_decay` | −8.086 | −0.304 | **× 6×10⁷** (channel switched on) |
| `k_ttca_deg` | −1.987 | −6.361 | ÷ 2×10⁴ |
| `k_fur_fft` | −4.257 | −0.576 | **× 4.8×10³** |
| `k_thiol_decay` | −3.045 | −5.976 | ÷ 850 |
| `k_thiolate_loss` | −2.056 | −1.051 | × 10.1 |
| `Ea_decay_carbonyl_sink` | 150.2 | **249.92** (250 ceiling) | traded onto its bound |
| `Ea_decay_thiol_sink` | 248.0 (250 ceiling) | 216.1 | traded off its bound |

These are topology switches, not refinements, and they were bought with a **10%
reduction in χ²ᵣ on a 10-degree-of-freedom objective explored from two starts.**
The wave's own diagnosis §9 records that B2.2's two starts gave 16.59 / 8.99 and
B2.3's gave 10.38 / 8.18 — a spread between starts of the same order as the
difference between waves. The exam is measuring which basin was found.

**No single revert restores the exam.** Reverting one constant at a time to its
B2.2 value, at B2.3's parameters, buffer-completed:

| revert | paired median | sulfur-only median |
|---|---:|---:|
| *(none — shipped)* | 50.13× | 130.6× |
| `k_thiolate_loss` | **42.17×** | **57.6×** |
| `k_pent_tdp` | 48.68× | 64.9× |
| `k_tdp_fur` | 52.62× | 68.6× |
| `decay_carbonyl_sink` | 50.13× | 100.7× |
| all others | 50.13× | 120–156× |

The best single revert recovers 0.08 of the 0.65 lost decades. The regression is
a **joint** re-routing, which is what a basin change looks like and what a single
bad constant does not.

---

## 4. The pH-response ladder — the row-level signature of the defect

Hofmann & Schieberle 1998 measured ribose+cysteine at 145 °C / 20 min in 0.5 M
phosphate at **pH 3, 5 and 7**. The pH-5 rung is a declared FIT anchor; pH 3 and
pH 7 are hold-out. Sweeping the declared label pH at frozen parameters:

| label pH | in-situ end | FFT B2.2 | FFT B2.3 | FFT measured | MFT B2.2 | MFT B2.3 | MFT measured |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3.0 | 3.38 | 333 | 1660 | **229** | 25.9 | 13.1 | **553** |
| 4.0 | 4.34 | 673 | 1823 | — | 378 | 308 | — |
| 4.5 | 4.68 | 500 | 1109 | — | 486 | 506 | — |
| **5.0** | 4.94 | 315 | 295 | **121** *(FIT)* | 427 | 340 | **198** *(FIT)* |
| 5.5 | 5.32 | 102 | **2.38** | — | 192 | **7.46** | — |
| 6.0 | 5.80 | 4.24 | 0.077 | — | 7.80 | 0.835 | — |
| 7.0 | 6.78 | 2.68 | **0.024** | **12** | 0.178 | **0.048** | **25** |

**pH3 → pH7, in decades:**

| | measured | B2.2 | B2.3 |
|---|---:|---:|---:|
| FFT | **+1.28** | +2.09 | **+4.84** |
| MFT | **+1.34** | +2.16 | +2.43 |

Three defects are visible at once, and only the third is new in B2.3:

1. **The shape is wrong in both waves.** The measured ladder is *monotone
   falling*; the model's is a **hump peaking near label pH 4–4.5**. On MFT the
   sign is inverted at the acid end — measured MFT is at its maximum at pH 3
   (553 ppb) where the model puts its minimum (13 ppb).
2. **The single fit anchor sits on the edge of a cliff.** Between label pH 5.0
   and 5.5 — in-situ 4.94 → 5.32, less than half a pH unit — B2.3's FFT falls
   2.1 decades. The fit is anchored at exactly one point on a curve whose shape
   nothing constrains, and it is anchored 0.4 in-situ units below a cliff.
3. **B2.3 made the alkaline limb 2.75 decades steeper than B2.2 on FFT, against
   a measured 1.28 decades total.** That is the regression, in one number.

The collapse is not arbitrary: the buffer-completed pH-7 pots end at in-situ
**6.59–6.70**, past the top rung of Kumazawa's retention ladder (pH 6.4, where
only 11% of spiked FFT survives 121 °C / 10 min), and are then run at 145 °C for
20 min. The model is extrapolating a *measured destruction* law past its last
rung and in temperature and in time, with nothing on the formation side pushing
back. **The tension is Kumazawa-vs-Hofmann and it is now sharp**: Kumazawa says
FFT is 89% destroyed at pH 6.4 in ten minutes at 121 °C; Hofmann says 12 ppb of
FFT is present after twenty minutes at 145 °C at pH 7. This is a different
conflict from the Hofmann-vs-Zhou MFT sign conflict the hold-out report already
records, and it should be named separately.

---

## 5. The panel did not stay level — its scorer cannot see what moved

**Four of the exam's 23 answered points are the same measurements as four
hold-out panel rows.** Not analogous rows — the same numbers:

| panel row | panel B2.3 fold | exam B2.3 fold | |
|---|---:|---:|---|
| `hofmann_ribose_pH3_FFT` | 7.23× | 7.25× | same measurement |
| `hofmann_ribose_pH3_MFT` | 0.0237× (= 42.2× inverted) | 42.23× | same measurement |
| `hofmann_ribose_pH7_FFT` | 0.00200× (= 499× inverted) | 498.99× | same measurement |
| `hofmann_ribose_pH7_MFT` | 0.00193× (= 517× inverted) | 517.27× | same measurement |

The panel registers exactly the collapse the exam registers. Its **gating count**
does not, because gating is a censored statistic:

| Hofmann pH axis | B2.1 | B2.2 | B2.3 |
|---|---|---|---|
| `hofmann_ribose_pH3_MFT` | 0.059× F | 0.047× F | 0.024× F |
| `hofmann_ribose_pH3_FFT` | 24.8× F | 1.47× **P** | 7.23× F |
| `hofmann_ribose_pH7_MFT` | 1.13× **P** | 0.0054× F | 0.0019× F |
| `hofmann_ribose_pH7_FFT` | 0.379× **P** | 0.208× F | 0.0020× F |
| **block score** | **2/4** | **1/4** | **0/4** |

`hofmann_ribose_pH7_FFT` fell **2.02 decades** and cost the gating count nothing,
because it was already outside the band. Across the whole panel:

* **22 of 34 gating rows were already failing in B2.2.** Of those, 11 got worse
  and 9 got better in B2.3: **5.44 decades lost, 4.02 gained, net +1.42 decades
  of degradation entirely invisible to the pass/fail count.**
* The panel's own *continuous* statistics moved the same way the exam did, just
  more gently: median |log₁₀ fold| **0.635 → 0.674** (4.31× → 4.72×), geometric
  mean fold **10.5× → 12.3×**.
* The worst invisible degradations: `hofmann_ribose_pH7_FFT` (+2.02 dec),
  `cerny2003_intact_skeleton_share_95C` (+1.72), `hofmann_ribose_pH7_MFT`
  (+0.45), `zhou_120C_dimer_share_pH6` (+0.42), `hofmann_ribose_pH3_MFT` (+0.30).

**So the scorecards agree.** One is a rank statistic over a lumpy 23-point pool
(amplifies), the other is a censored count (suppresses). Neither is measuring
accuracy. The diagnosis's contradiction 1 — "the two scorecards disagree about
this wave, and nothing here reconciles them" — can be retired and replaced with
"the exam's median and the panel's gating count are both unstable summaries of
the same degradation, which is concentrated on the Hofmann pH-3/pH-7 axis."

*(Bookkeeping note: the hold-out JSON flags **34** rows `gating_in_b2_2` against
a headline denominator of **33**. `b2_2_pass` is 12 either way, so nothing here
turns on it, but it should be reconciled.)*

---

## 6. The five pH-vs-source disagreements account for exactly zero

Upper bound, not an estimate. The five bundles §7 of the diagnosis lists —
Bi 2020 ×2, Liu 2023, Li 2026, Steinhagen 2021 — carry **nine** exam points
between them, and **all nine are DECLINED by the core in both columns**. They
carry no fold error, enter no median, and are not in the paired subset. Four of
the five are `matrix_path_lipid`, which the core declines for having no
lipid-oxidation path at all; the fifth (Steinhagen HMF) is declined too.

They remain a real data defect — Steinhagen's `ph: 4.36` being an *after*-
autoclaving reading on a scale where every other bundle is an initial bench
reading is the serious one, and it would bite the moment that bundle became
answerable — but **they have zero leverage on this regression.**

---

## 7. A reproducibility defect in the pre-registration's own evidence

`kinetic_core_b2_3_prereg.md` §0 discloses a pre-freeze probe,
`scratch/b23_encoding_probe.py`, as the evidence for `AMINE_FATE_BASIS` — the
wave's self-declared "strongest single assumption" — and prints a three-row
table comparing three encodings.

**That probe does not run on the shipped tree.** It raises
`KeyError: 'AMN'` on its first case. `AMN` is not a species anywhere in
`src/kinetic_core/`: the shipped encoding does not *destroy* an amine pool, it
**has no amine pool**. Consequently two of the probe's three "encodings"
(`carry_carboxyl_only` and `carry_both`) reduce to the same unpatched code path,
and D1's own attempt to re-probe the axis at frozen B2.3 parameters was silently
a no-op — `enc_carry_both` came back bit-identical to the shipped cell on all
eighteen sulfur rows.

The most likely history is benign: an `AMN` species existed in the working tree
when the probe ran and was removed once the choice was made. That is an honest
sequence with a reproducibility defect, not a fabricated table. But the
consequence stands: **the amine-fate axis cannot currently be re-probed**, the
only encoding axis live in the tree is `CBX` on/off, and that axis is worth
nothing on the exam (Hofmann median 12.21× carried vs 12.02× not carried).

---

## 8. VERDICT

**The defect is not one defect; it is two, and neither is the one the wave was
sent to fix.**

> **D1-A — THE OBJECTIVE PRICES A pH ENDPOINT AT NINE TIMES A LEVEL MISS AND HAS
> NO USABLE CONSTRAINT ON FORMATION-VS-pH.** Three `zhou_final_pH_*` rows at
> σ = 0.25 *pH units* sit in the same sum of squares as 55 rows at
> σ_log = 0.2–0.6 *log-fold*. One of them, `zhou_final_pH_from_pH8`, carries 44%
> of the entire refit. Meanwhile the only two rows that measure how thiol
> *formation* responds to pH are declared unfittable in the fit script's own
> note, and every Hofmann anchor in the fit is pH 5.0. The optimiser was
> therefore free to buy two pH endpoints and one Kang pot by re-routing the
> network through 3–8-decade moves in `k_fur_fft`, `k_osone_decay`,
> `k_dimer_decay` and `k_ttca_deg`, and to pay for it on an axis nothing scored.
> **Evidence:** §3. **Predicted signature if fixed:** the pH3→pH7 FFT slope falls
> from +4.84 decades back toward the measured +1.28; the four Hofmann pH-3/pH-7
> rows come off 318–517× and the panel's Hofmann block returns to ≥1/4; the
> exam's geometric-mean fold falls from 47.7× toward B2.2's 26.8×.

> **D1-B — BOTH SCORECARDS' HEADLINE STATISTICS ARE UNSTABLE, AND THEY SHARE
> ROWS.** The exam's paired median is a rank statistic over 23 heterogeneous
> points whose reported value in both waves is a row that did not move; four
> crossings swing it 4.7× while the geometric mean moves 1.78×. The panel's
> gating count is censored and hid 1.42 net decades. And four points appear in
> *both*, so the two are not independent evidence. **Evidence:** §2, §5.
> **Predicted signature if fixed:** the two scorecards stop appearing to
> disagree, because both will be reporting the same continuous quantity.

**On the mandate's sharpest question — is a true 50× worth more than a false
10×?** The premise is half right and half wrong, and the distinction matters:

* **B2.2's 10.65× was not false, and it was not propped up by the charge leak.**
  The `S_only` cell settles this: with the leak repaired and B2.2's parameters
  kept, the exam median is 10.63× and the in-band count is unchanged at 5/23.
* **B2.3's 50.13× is not a "truer" measurement of the same model either.** It is
  the same pool re-sorted by a genuinely worse parameter vector, and the honest
  size of the degradation is the geometric mean: **26.8× → 47.7×, a factor of
  1.78**, not 4.7.
* **The finding that survives both corrections is the uncomfortable one:** the
  refit improved in-sample χ²ᵣ (1.80 → 1.635) and made out-of-sample accuracy
  worse by every statistic on both scorecards, on a 58-row / 48-parameter /
  10-dof objective explored from two starts. That is textbook overfitting, and
  the charge-conservation fix — which is sound, free, and vindicated here — is
  not responsible for it.

### What a B2.4 should change

Not a stoichiometry basis, not a centre delta, not an engine input. In order:

1. **A DECLARED WEIGHTING.** Make the pH-unit-vs-log-fold exchange rate an
   explicit, defended declaration rather than an emergent consequence of two
   independently chosen σ. Then run the fit at 2–3 declared weightings and
   **report the exam under each**, rather than picking one and reporting a
   single number.
2. **MULTI-START DISCIPLINE, OR FEWER FREE PARAMETERS.** Two starts on a 48-dim,
   10-dof objective does not license calling the reported point *the* optimum,
   and the exam is measuring basin choice. Either report an ensemble (N starts,
   exam median with spread) or freeze the constants the corpus does not identify
   — 37 of 48 are unidentified even at B2.3.
3. **A SCORER CONDITION, ON BOTH SCORECARDS.** The exam should lead with the
   geometric-mean fold and per-family medians, keeping the paired median as a
   secondary. The hold-out panel should print median |log₁₀ fold| beside the
   gating count so a 100× degradation on an already-failing row is visible. And
   the four shared rows must be declared as shared.
4. **RETIRE OR REPAIR `scratch/b23_encoding_probe.py`**, and restate
   `AMINE_FATE_BASIS` as what the code does — the amine is not represented —
   rather than as a choice among three encodings the tree cannot express.
5. **NAME THE KUMAZAWA-VS-HOFMANN CONFLICT** as a corpus contradiction in its
   own right, separate from the Hofmann-vs-Zhou MFT sign conflict already
   recorded. No parameter choice satisfies both at pH ≥ 6.5.

### Drafted pre-registration text for B2.4

> **B2.4 — PRE-REGISTRATION (draft, D1)**
>
> **What this wave changes.** No stoichiometry, no species, no benchmark, no
> engine input. Two things only: (i) `generate_kinetic_core_b2_*_fit.py` gains a
> single declared scalar `PH_ENDPOINT_WEIGHT` with a written basis, applied to
> the three `zhou_final_pH_*` rows, and the fit is run at three declared values
> of it; (ii) the fit is run from **six** seeded starts at each weighting rather
> than two, and the ensemble — not the best member — is what is reported.
>
> **W-1.** At the shipped weighting, six starts produce a spread in final cost of
> **more than 2.0** (B2.2 saw 16.59/8.99, B2.3 10.38/8.18), confirming that the
> two-start protocol was selecting a basin rather than finding an optimum.
> *Falsified if the six starts converge to within 0.5.*
>
> **W-2.** Across the ensemble at the shipped weighting, the exam's
> **geometric-mean fold** spans a range of **at least 1.5×**. If it does, the
> single number "50.13×" was never a property of the model.
>
> **W-3.** Down-weighting the pH endpoints by 4× moves `k_fur_fft` back below
> −2.0 (B2.3: −0.576; B2.2: −4.257) and the pH3→pH7 FFT slope below **+3.0
> decades** (B2.3: +4.84; measured: +1.28).
>
> **W-4.** At every weighting tried, the **charge-conservation fix stays free**:
> the exam's geometric-mean fold at B2.2's parameters with B2.3's stoichiometry
> stays within 10% of 26.8×. This wave predicts its predecessor's conservation
> fix is not the problem and says so in advance.
>
> **W-5.** The Hofmann pH-3/pH-7 hold-out block does **not** reach 3/4 at any
> weighting. The shape defect in §4 — a hump against a monotone ladder — is not a
> weighting problem and re-weighting will not fix it. *This is the wave's
> prediction that its own remedy is insufficient.*
>
> **W-6.** `fed_c2c3_MFT_pH3` remains the largest single residual in the fit at
> every weighting, because the route it scores still carries no pH factor.
>
> **W-7.** The exam and the hold-out panel, once both report median |log₁₀ fold|,
> agree in **sign** on every weighting tried. If they ever disagree in sign, the
> shared-row overlap is masking something and the wave says so.
>
> **What would make B2.4 a failure.** W-1 falsified (the objective was not
> multi-modal and D1's diagnosis was wrong); or the ensemble's exam spread being
> narrow while the between-wave difference stays large, which would mean the
> regression is a real physics change after all and D1 mis-attributed it.
