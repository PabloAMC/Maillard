# Directional (ordinal) accuracy of the Maillard model

**Wave S2, 2026-08-27.** Panel: [`directional_claims_panel.yml`](directional_claims_panel.yml) (52 claims, 14 sources).
Runner: scratchpad only — deliberately not committed, so no future refactor can quietly
start optimising against it. Nothing in the model was tuned to produce these numbers.

> **Provenance.** Scored against `HEAD = c1a12d2` with 5 uncommitted files in the tree (a
> concurrent Wave S agent was editing `src/recommend.py` during this run). The panel was
> scored twice, before and after that agent's commit, and all 52 rows were identical — the
> headline is not an artifact of a moving tree. The panel file is `.yml`, **not** `.yaml`,
> because `.gitignore:117` ignores `*.yaml` wholesale (a Cantera-output rule, with
> exceptions only under `data/`); a `.yaml` panel would have been silently untracked.

---

## 1. The question

Every existing validation artifact in this repo measures **absolute** accuracy, and the
honest answer there is bad: 0/8 on true extrapolation for the matrix lane, 36x
over-prediction on lipid aldehydes, order-of-magnitude scatter elsewhere. But nobody uses
a model like this to predict a ppb number. They use it **ordinally**:

> "Should I switch from glucose to ribose?" · "Which of these two isolates will taste
> beanier?" · "Which way do I move the barrel temperature?"

Those are sign and ranking questions. Absolute error is close to irrelevant to them. This
document is the first measurement of whether the model gets **directions** right.

---

## 2. Headline

> **This is the Wave S2 measurement, not the current one.** Two later waves moved rows. The
> shipped tree scores **21/29**, not the 18/29 below. Jump to
> [CURRENT STANDING](#current-standing--the-panel-as-it-scores-on-the-shipped-tree) at the end
> of this file for the numbers that are live, and read this section as the baseline the delta
> is measured from. The pre-fix numbers are kept verbatim because deleting them would delete
> the delta.

| Bucket | Agree | Evaluable | Rate | Unevaluable |
|---|---|---|---|---|
| **Strictly independent (headline)** | **18** | **29** | **62%** | 8 |
| + system-overlap rows | 6 | 6 | 100% | 0 |
| **Screening total (ind + overlap)** | **24** | **35** | **69%** | 8 |
| Fit-adjacent (**excluded from headline**) | 9 | 9 | **100%** | 0 |

Read the last two rows together. The model agrees with **100%** of the claims whose
benchmarks it was fitted against, **100%** of the claims that share a chemical system with
a fitted row, and **62%** of the claims that are genuinely out of sample. That gradient is
the whole story: the model is a good interpolator of what it was shown and a mediocre
extrapolator, and the ordinal score degrades exactly where you would predict.

62% is above chance but not by as much as it looks. Most panel rows are binary (A>B vs
A<B), so a coin flip scores ~50%. **The model's directional edge over a coin is real but
modest — and it is not spread evenly.**

## 3. Where the failures live

| Category | Agree / evaluable | Verdict |
|---|---|---|
| sugar_identity | **8 / 8** | Trustworthy |
| temperature | **6 / 8** | Mostly trustworthy; fails on non-monotonic |
| time | 2 / 2 | Trustworthy (sign only — see §6) |
| additive_cysteine | 3 / 4 | Trustworthy |
| lipid_lane | 2 / 2 | Trustworthy |
| matrix_identity | 1 / 1 | Thin (one row) |
| **ph** | **2 / 7** | **Broken — worse than a coin** |
| **moisture_aw** | **0 / 3** | **Broken — no usable signal** |

Two knobs carry almost all the failure. Drop pH and aw and the headline is **16/19 (84%)**.
Keep only pH and aw and it is **2/10 (20%)**, which is *worse than chance* — the model is
anti-correlated with the literature on those axes. That is more useful than a flat 62%:
it says the model is a usable screening tool over sugar identity, precursor loading,
temperature and time, and an actively misleading one over pH and moisture.

### 3.1 pH — the model's pH physics does not reach the prediction path

Two independent sources, measured directly, say pyrazines RISE with pH:

* Laemont & Barringer 2023 (roasted sunflower, Table 2): dimethylpyrazine isomers
  26.6 → 37.4 → 68.2 ppb across pH 4 → 7 → 9. **Model: 99.5 → 16.8 → 14.9** (PH-06).
* Cerny 2015, quoting Farmer & Mottram: "Higher pH produced more thiazoles"; "At pH 4.5,
  sulfur compounds dominated the volatile fraction, whereas with pH 6.5 nitrogen compounds
  prevailed." **Model: 82.4 → 41.6 → 25.1** across pH 4.5 → 5.6 → 6.5 (PH-04).

The model moves pyrazines the wrong way, monotonically, in two different systems. Three
things in the code explain it, all verifiable by inspection:

1. **`ReactionConditions.get_ph_multiplier` is never called on the prediction path.** It is
   the function that encodes the documented physics — "2,3-enolisation (pyrazines):
   Favored at alkaline pH", a 0.2→8.2 sigmoid, and a 4.9x acidic boost for the furan/thiol
   track. `grep -rn get_ph_multiplier src/` finds it only in `kinetics.py`,
   `pathway_ranker.py` and `cantera_export.py`. `benchmark_validation._run_benchmark_recommendation`
   calls `conditions.get_rate_constant(...)`, which applies `_ionization_correction` and
   `_water_activity_correction` and **not** `get_ph_multiplier`. The pH physics is written,
   unit-tested (`tests/unit/test_conditions.py` asserts the correct signs on it), and dead.
2. **The pyrazine branch of `_ionization_correction` is unreachable.** It keys on the
   substring `"pyrazine"` in the reaction-family name (`conditions.py:310`). The families
   the SMIRKS engine actually emits for a glycine/glucose system are
   `Aminoketone_Condensation`, `Enolisation_2_3`, `Enolisation_2_3_Amadori`,
   `Retro_Aldol_Fragmentation`, … — **none contains "pyrazine"**, so every one of them
   returns exactly 1.0 at every pH.
3. **The only live pH term is the amine-deprotonation factor** (pKa 8.0) applied to
   `amadori` / `schiff` / `strecker`, which runs 1e-4 at pH 4 to 0.91 at pH 9. It boosts
   the shared condensation trunk ~9000x across that range; the pyrazine and furanone
   branches then lose share in the projection budget. What survives at the output is a
   redistribution artifact, not chemistry.

The same mechanism inverts furfural: the code intends acidic-favoured (a 4.9x boost at
pH 4), and the output rises 403 → 917 across pH 4 → 7 (PH-05). Notably, when the panel
uses a *measured* furfural-vs-pH series instead of a textbook inference, the literature
says furfural is **flat** — "Furfural was also not affected by a change in pH" — and the
model still misses, by rising 21% (PH-07). PH-05 and PH-07 are both scored, and both miss;
PH-07 is the better-sourced of the two.

### 3.2 Water activity — implemented, correctly shaped, and inert

`_water_activity_correction` does implement the Labuza curve and does peak at aw 0.65.
Measured directly:

```
aw     Amadori/Strecker   Enolisation_1_2   Enolisation_2_3   Schiff   Aminoketone   Retro_Aldol
0.25       0.200              1.000             1.000          1.000      1.000         1.000
0.65       1.000              1.000             1.000          1.000      1.000         1.000
0.98       0.623              1.000             1.000          1.000      1.000         1.000
```

The correction is applied to **2 of the 7 families** the engine emits, and neither of them
is on the furfural/HMF track. By the time the signal reaches an observable it is a 4.9%
wiggle: HMF over aw 0.25 → 0.65 → 0.95 comes out 570 → 598 → 594. The *shape* is right —
it does peak at 0.65 — but the amplitude is below the panel's 5% no-separation threshold,
so all three aw rows score as "no signal" (AW-01, AW-03) or invert (AW-02). A formulator
asking "should I run this drier?" gets nothing back.

This is a better diagnosis than "aw is broken": the physics is present and correctly
shaped, and the fix is a routing question (which families the correction is applied to),
not a modelling one.

### 3.3 The two genuine chemistry misses

* **TEMP-01, acrylamide vs temperature.** Ma 2024 measures 130 > 150 > 170 °C in extruded
  PBMA ("the highest amount observed at ... 130 °C"); the model gives 0.63 → 4.4 → 31.6,
  monotonically rising, a 50x swing in the wrong direction. Acrylamide degrades and
  re-adds above ~160 °C and the model has no destruction term at all — confirmed by
  ACR-02, where a review states the acrylamide peak sits at 160–180 °C and the model runs
  straight through it (0.70 → 38.7 → 1937 over 130 → 170 → 210 °C). **The model can only
  make acrylamide go up with temperature.** For a tool whose selling point includes
  process-safety screening, this is the most consequential single failure in the panel.
* **CYS-02, cysteine suppressing the browning track.** Model moves HMF by −3.9%, below the
  separation threshold. Directionally right, magnitude negligible.

## 4. Where it is genuinely good

**Sugar identity is 8/8, and 6 of those 8 are strictly independent rows sourced from
papers with no relationship to any benchmark in this repo.** Ni 2021 (camellia hydrolysate,
110 °C) independently reproduces the pentose > hexose ordering for furfural, MFT and FFT
that Hofmann & Schieberle reported, and the model gets all three. It also gets the
class-total pyrazine ordering ribose > glucose > fructose (SUG-10), the fructose > glucose
HMF ordering from a baking system (SUG-12), and — the one the panel expected it to lose —
the *inverted* glucose > fructose ordering for furfural (SUG-13).

Two others worth naming:

* **TEMP-04** is the single most on-target row in the panel: a cornmeal / pea-protein-isolate
  extrudate, a barrel-temperature contrast (145 → 160 °C max zone), and a named model
  observable. The source says pyrazines "nearly doubled"; the model says 0.74 → 2.80. Right
  sign, roughly 2x too much magnitude.
* **ACR-01** gets a genuine non-monotonic right: acrylamide peaking at intermediate
  asparagine-source loading because the reducing-sugar co-substrate is diluted out. The
  model reproduces the peak (0.28 → 0.77 → 0.54). So the machinery *can* express
  non-monotonic behaviour — it is specifically the temperature and aw axes where it cannot.

## 5. Full panel

`ind` = strictly independent (headline) · `overlap` = shares a system with a fitted row
(ratio the fit cancels out of) · `FIT` = declared per-row fit target, **excluded from the
headline**.

| Claim | Category | Observable | Literature | Model | Outcome | Fit | Predicted values (ppb) |
|---|---|---|---|---|---|---|---|
| SUG-01 | sugar_identity | MFT | `A>B` | `A>B` | agree | FIT | cysteine + D-ribose = 242.4; cysteine + D-glucose = 74.27 |
| SUG-02 | sugar_identity | FFT | `A>B` | `A>B` | agree | FIT | cysteine + D-ribose = 218; cysteine + D-glucose = 139.4 |
| SUG-03 | sugar_identity | FFT | `A>B` | `A>B` | agree | ind | FFT in cysteine + glucose = 139.4; MFT in cysteine + glucose = 74.27 |
| SUG-04 | sugar_identity | MFT | `A>B` | `A>B` | agree | FIT | MFT in cysteine + ribose = 242.4; FFT in cysteine + ribose = 218 |
| SUG-05 | sugar_identity | Furfural | `A>B` | `A>B` | agree | overlap | cysteine + D-ribose = 511.9; cysteine + D-glucose = 325.4 |
| SUG-06 | sugar_identity | MFT | `A>B` | `-` | n/a | ind | - |
| PH-01 | ph | MFT | `A>B` | `A>B` | agree | overlap | pH 4.5 = 271.2; pH 6.5 = 40.23 |
| PH-02 | ph | bis(2-methyl-3-furyl) disulfide | `A>B` | `A>B` | agree | overlap | pH 4.5 = 219.3; pH 6.5 = 16.68 |
| PH-03 | ph | FFT | `decreasing` | `increasing` | **MISS** | ind | pH 4.0 = 0.0009; pH 5.5 = 0.0396; pH 7.0 = 0.0379 |
| PH-04 | ph | 2,5-Dimethylpyrazine | `increasing` | `decreasing` | **MISS** | ind | pH 4.5 = 82.42; pH 5.6 = 41.64; pH 6.5 = 25.14 |
| PH-05 | ph | Furfural | `decreasing` | `increasing` | **MISS** | ind | pH 4.0 = 402.9; pH 5.5 = 741.1; pH 7.0 = 917.5 |
| TEMP-01 | temperature | Acrylamide | `decreasing` | `increasing` | **MISS** | ind | 130 C = 0.6341; 150 C = 4.405; 170 C = 31.55 |
| TEMP-02 | temperature | Acrylamide | `A>B` | `A>B` | agree | ind | extruded, 130 C 20 min = 0.6341; unheated control, 25 C = 0.0004 |
| TEMP-03 | temperature | HMF | `increasing` | `increasing` | agree | ind | 130 C = 184.5; 150 C = 997.7; 170 C = 4625 |
| AW-01 | moisture_aw | HMF | `decreasing` | `flat` | **MISS** | ind | aw 0.3 = 580.3; aw 0.6 = 597.2; aw 0.98 = 593 |
| AW-02 | moisture_aw | Acrylamide | `peak` | `trough` | **MISS** | ind | aw 0.3 (low moisture) = 8.042; aw 0.6 (mid moisture) = 3.904; aw 0.9 (high moisture) = 4.359 |
| CYS-01 | additive_cysteine | FFT | `A>B` | `A>B` | agree | ind | glycine + glucose + cysteine = 33.11; glycine + glucose (no cysteine) = 0 |
| CYS-02 | additive_cysteine | HMF | `A>B` | `flat` | **MISS** | ind | glycine + glucose (no cysteine) = 593; glycine + glucose + cysteine = 570 |
| SCOPE-01 | scope | CEL | `A>B` | `-` | n/a | ind | - |
| SCOPE-02 | scope | Methional | `A>B` | `-` | n/a | ind | - |
| MAT-01 | matrix_identity | Hexanal | `A>B` | `-` | n/a | ind | - |
| MAT-02 | matrix_identity | Hexanal | `A>B` | `A>B` | agree | FIT | soy isolate, 40 C = 1297; pea isolate, 40 C = 889.6 |
| MAT-03 | matrix_identity | Hexanal | `A>B` | `A>B` | agree | overlap | hexanal, pea isolate 40 C = 889.6; 1-hexanol, pea isolate 40 C = 80.1 |
| PROC-01 | process_heating | Hexanal | `A>B` | `A>B` | agree | FIT | processed, 140 C = 771.1; unprocessed control, 25 C = 1.839 |
| PROC-02 | process_heating | 2-Pentylfuran | `A>B` | `A>B` | agree | FIT | processed, 140 C = 160.7; unprocessed control, 25 C = 1.043 |
| PROC-03 | process_heating | Nonanal | `A>B` | `A>B` | agree | FIT | processed, 140 C = 10.41; unprocessed control, 25 C = 0.1371 |
| PROC-04 | process_heating | 2,5-Dimethylpyrazine | `flat` | `-` | n/a | ind | - |
| PROC-05 | ranking | Hexanal | `A>B` | `A>B` | agree | FIT | hexanal = 771.1; 2-pentylfuran = 160.7 |
| TIME-01 | time | HMF | `increasing` | `increasing` | agree | ind | 5 min = 100.4; 30 min = 593; 120 min = 2334 |
| LIP-01 | lipid_lane | Hexanal | `A>B` | `A>B` | agree | FIT | hexanal, pea isolate (linoleic 50%, oleic 22%) = 800.6; nonanal, pea isolate = 10.81 |
| LIP-02 | lipid_lane | Hexanal | `increasing` | `increasing` | agree | overlap | 0.1 min = 800.6; 10 min = 8.006e+04; 60 min = 4.804e+05 |
| LIP-03 | lipid_lane | Hexanal | `A>B` | `-` | n/a | ind | - |
| SUG-07 | sugar_identity | Furfural | `A>B` | `A>B` | agree | ind | cysteine + D-ribose, 110 C 90 min = 184.2; cysteine + D-glucose, 110 C 90 min = 89.3 |
| SUG-08 | sugar_identity | MFT | `A>B` | `A>B` | agree | ind | cysteine + D-ribose, 110 C 90 min = 5.286; cysteine + D-glucose, 110 C 90 min = 0.8398 |
| SUG-09 | sugar_identity | FFT | `A>B` | `A>B` | agree | ind | cysteine + D-ribose, 110 C 90 min = 4.599; cysteine + D-glucose, 110 C 90 min = 2.63 |
| SUG-10 | sugar_identity | 2,5-Dimethylpyrazine | `decreasing` | `decreasing` | agree | ind | D-ribose = 3.127; D-glucose = 1.748; D-fructose = 0.4613 |
| SUG-11 | sugar_identity | FFT | `A<B` | `-` | n/a | ind | - |
| SUG-12 | sugar_identity | HMF | `A>B` | `A>B` | agree | ind | glycine + D-fructose, 150 C 35 min = 2433; glycine + D-glucose, 150 C 35 min = 1524 |
| SUG-13 | sugar_identity | Furfural | `A>B` | `A>B` | agree | ind | glycine + D-glucose, 165 C 8 min = 861; glycine + D-fructose, 165 C 8 min = 177.3 |
| PH-06 | ph | 2,5-Dimethylpyrazine | `increasing` | `decreasing` | **MISS** | ind | pH 4 = 99.46; pH 7 = 16.76; pH 9 = 14.94 |
| PH-07 | ph | Furfural | `flat` | `increasing` | **MISS** | ind | pH 4 = 751.8; pH 7 = 901.9; pH 9 = 908.1 |
| TEMP-04 | temperature | 2,5-Dimethylpyrazine | `A>B` | `A>B` | agree | ind | T4 profile, 160 C max = 2.797; T1 profile, 145 C max = 0.7392 |
| TEMP-05 | temperature | HMF | `increasing` | `increasing` | agree | ind | 150 C = 1096; 170 C = 4986; 190 C = 1.945e+04 |
| TEMP-06 | temperature | Furfural | `A>B` | `A>B` | agree | ind | roasted, 165 C 8 min = 901.9; raw, 25 C = 0.4482 |
| AW-03 | moisture_aw | HMF | `peak` | `flat` | **MISS** | ind | aw 0.25 = 570.5; aw 0.65 = 598.2; aw 0.95 = 593.7 |
| CYS-03 | additive_cysteine | MFT | `A>B` | `A>B` | agree | ind | hydrolysate + xylose + cysteine, 100 C 60 min = 0.1408; hydrolysate + xylose, no cysteine = 0 |
| CYS-04 | additive_cysteine | 2,5-Dimethylpyrazine | `A>B` | `A>B` | agree | ind | hydrolysate + xylose (no cysteine), 120 C = 2.215; hydrolysate + xylose + cysteine, 120 C = 1.742 |
| LIP-04 | lipid_lane | Hexanal | `A>B` | `A>B` | agree | overlap | hexanal, soy isolate (linoleic 53%) = 6.99e+05; nonanal, soy isolate (oleic 23%) = 1.932e+04 |
| ACR-01 | temperature | Acrylamide | `peak` | `peak` | agree | ind | low protein (asn 3 mM, glucose 10 mM) = 0.2847; mid protein (asn 10 mM, glucose 7 mM) = 0.7719; high protei... |
| ACR-02 | temperature | Acrylamide | `peak` | `increasing` | **MISS** | ind | 130 C = 0.7009; 170 C = 38.73; 210 C = 1937 |
| TIME-02 | time | HMF | `increasing` | `increasing` | agree | ind | 5 min = 222.9; 25 min = 1096; 35 min = 1524 |
| SCOPE-03 | scope | 2-Pentylfuran | `A>B` | `-` | n/a | ind | - |

## 6. The fit-adjacent rows, reported separately

Nine claims are excluded from the headline because a declared per-row fit target is
literally one side of the comparison. All nine agree. This is **not** evidence about the
model:

| Claim | Why excluded | Fit record |
|---|---|---|
| SUG-01, SUG-02, SUG-04 | A-side (or both sides) is `cys_ribose_140C_Hofmann1998` | `sulfur_barrier_refit_hofmann.json` — 4 free barriers / 2 rows, `per_row_recovery` |
| MAT-02 | Both sides are the fitted pea/soy 40 °C rows | `matrix_observability_refit_pratap_singh.json` — 1 shared scale / 2 rows, `per_row_recovery` |
| PROC-01, PROC-02, PROC-03, PROC-05, LIP-01 | A-side is `pea_isolate_uht_140C_Trikusuma2019`; its declared targets ARE 782 / 163 / 24 ppb | `projection_constant_refit.json` |

PROC-01 makes the point concretely: the model returns hexanal = 771 ppb against a
benchmark target of 782 ppb. That is recovery to 1.4%, and it tells you nothing about
whether the model would have got the direction right had it not been shown the answer.

The `overlap` tier is a softer version of the same worry and is reported separately rather
than excluded, because these are ratios the fit's free parameters cancel out of — a single
shared scale fitted at pH 5 cannot determine a pH-4.5-vs-6.5 ratio. All 6 agree. If you
disbelieve that argument, use the strict 18/29; if you accept it, use 24/35. Both are in
the table and neither is the number I would quote alone.

## 7. What the model provably cannot say (8 unevaluable rows)

These are not gaps in the literature search. Each was probed against the running model and
the limitation measured.

| Claim | The limitation |
|---|---|
| SUG-06, SUG-11 | **The model cannot distinguish D-ribose from D-xylose.** Substituting one for the other gives identical predictions to 4+ significant figures on every observable. Ni 2021 measures xylose giving **2.3x** the FFT of ribose while ribose gives **1.8x** the MFT of xylose — a large, sign-crossing real difference. Any pentose recommendation the model makes is a coin flip on this axis. |
| SCOPE-01 | **No α-dicarbonyls, no CML, no CEL.** Ma 2024's central finding is that CML and CEL move in *opposite* directions with both temperature and moisture. Neither compound, nor glyoxal / methylglyoxal / 3-deoxyglucosone, appears in the model's output for any probed precursor set — despite `cml_cel_commercial_pbma_Foods2023.json` existing as a benchmark. |
| SCOPE-02 | **No methionine Strecker branch.** L-Methionine resolves as a precursor, but methionine + D-glucose at 140 °C / pH 6 / 30 min returns **zero** named volatiles: no methional, no dimethyl disulfide. Methional is one of the largest relative responses in Trikusuma's UHT table (0.55 → 3.10 µg/L). |
| SCOPE-03 | **2-Pentylfuran has no Maillard route in the model, only a lipid one.** Ni 2021 measures hexoses giving *more* 2-pentylfuran than pentoses — the opposite sign to the thiols. The free-precursor lane never emits 2-pentylfuran for any sugar, so this sign cannot be represented at all. |
| MAT-01 | **`MATRIX_BENCHMARK_PROFILES` has exactly two entries**, `pea_iso` and `soy_iso`. Brown rice — which in Pratap-Singh's Table 1 carries **14–20x** the hexanal of either legume isolate, by far the largest effect in the table — cannot be placed in the ranking. |
| PROC-04 | **The matrix lane and the free lane share no observable.** For a protein isolate the model emits only Hexanal, 2-Pentylfuran, 1-Hexanol and Nonanal — no pyrazine, no methional, nothing on the Maillard track. So it cannot answer "what happens to pyrazines in a pea beverage", which is the question a formulator of that product would actually ask. |
| LIP-03 | **No off-note can go DOWN with heat.** Trikusuma measures "beany" falling 3.3 → 1.5 while "oxidized/painty" rises 1.0 → 2.7 under the same UHT step. The model's four matrix markers are all oxidised-green descriptors and all rise with thermal load. It has no representation of a thermal step *destroying* an off-note — which is the entire mechanism behind the most common real-world pea-protein intervention. |

## 8. What this score does and does not license

**It licenses screening claims of a specific, narrow shape.** If the question is "which of
these two sugars, or which of these two temperatures, or which of these two precursor
loadings gives more of compound X", and X is one of the nine free-lane or four matrix-lane
observables, and pH and moisture are held FIXED between the two arms, then the model is
right about 8 times in 10 (16/19 on those axes). That is genuinely useful for triaging a
formulation grid down to a shortlist before anyone runs an experiment. It is the model's
best-supported use case and, until now, an unmeasured one.

**It does not license quantitative claims.** Nothing here contradicts the absolute-accuracy
findings; the panel deliberately never scores a magnitude. Where magnitudes are visible
they are frequently poor even when the sign is right — SUG-13 gets glucose > fructose for
furfural but by 4.9x where the measurement is 1.2x; TEMP-04 gets the right sign with about
twice the reported effect. **A correct sign here is not evidence for a correct number.**

**It does not license pH or moisture recommendations, at all.** 2/10 on those axes is
worse than a coin, and the failures are systematic rather than noisy: the model moves
pyrazines the wrong way with pH in every system tested, and has no usable aw response. A
user who moves pH on this model's advice is more likely to be misled than helped. This
should be a runtime guard, not a footnote.

**It does not license any claim about acrylamide and temperature.** The model is
monotonically increasing where two sources place a maximum in the process-relevant window.

**62% is not a quality certificate.** Most rows are binary, so chance is ~50%. The honest
one-line summary is: *the model has a real but modest directional edge, concentrated
almost entirely in sugar identity, precursor loading, temperature and time, and it is
anti-informative on pH and water activity.*

**Two caveats on the panel itself.**
*Selection.* The claims were chosen to fit what the model can express. Categories where it
emits nothing (α-dicarbonyls, Strecker aldehydes, brown rice) are recorded as unevaluable
rather than as failures — which is the honest bookkeeping, but it means the panel is
scoped to the model's own footprint and cannot be read as coverage of Maillard chemistry.
*Condition mapping.* Several literature systems are food matrices (cookies, extrudates,
roasted seeds) mapped onto a two-precursor aqueous model at an assumed pH and aw. Those
mappings are recorded per claim; where one is loose the row carries a `fit_note` saying so.
A miss on a loosely-mapped row is weaker evidence than a miss on TEMP-01 or PH-06.

## 9. Reproducing this

The runner lives in the session scratchpad, not the repo. This is deliberate: an in-repo
directional scorer becomes an optimisation target the moment someone runs it in CI, and the
value of this panel is that the model has never seen it. To re-run:

```bash
source /opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh && conda activate maillard
export PYTHONPATH=/Users/pabloantoniomorenocasares/Developer/Maillard
python <scratchpad>/run_directional_panel.py     # ~40 s, writes directional_results.json
```

The runner reads `docs/validation/directional_claims_panel.yml`, builds a synthetic
benchmark payload per condition, and calls `src.benchmark_validation.evaluate_benchmark_payload`
— the same entry point the benchmark suite uses. Every quote in the panel was verified
programmatically against the retrieved source text; the 14 quotes that are not contiguous
verbatim strings carry a `quote_form` field explaining exactly what was changed.

**If you re-run this after changing the model, report the strictly-independent number and
the per-category breakdown together.** A single headline that rose while pH stayed at 2/7
would be hiding the finding.

---

# Appendix — Mission 2: three hard retrievals

## 2a. Trikusuma — **CONFIRMED**

`data/benchmarks/pea_isolate_uht_140C_Trikusuma2019.json` claims hexanal 782 / 2-pentylfuran
163 / nonanal 24 ppb. All three are real.

**Source found:** Trikusuma, M., *Identification of Aroma Compounds in Pea Protein UHT
Beverages*, MS thesis, The Ohio State University, 2018 (advisor D.G. Peterson) — the direct
precursor of the paywalled Food Chemistry paper. OhioLINK ETD accession
`osu1531495328317918`, open access, no auth. Table 6, p.35, columns
`Control ± Stdev (µg/L) | Processed ± Stdev | Aged ± Stdev | Odor Threshold`:

```
 802   Cut grass              Hexanal          331.00 ± 81.27a    781.72 ± 58.59b    683.08 ± 58.09b     4.5
1003   Solvent-like, floral   2-Pentylfuran     59.39 ± 1.93      163.16 ± 15.06     197.03 ± 6.01       6
1103   Plastic, waxy          Nonanal            8.24 ± 0.44       23.98 ± 0.80b      22.44 ± 0.28       1
```

| Benchmark | Thesis | |
|---|---|---|
| hexanal 782 ppb ±8% | 781.72 ± 58.59 µg/L (7.5% RSD) | match |
| 2-pentylfuran 163 ppb ±9% | 163.16 ± 15.06 (9.2% RSD) | match |
| nonanal 24 ppb ±4% | 23.98 ± 0.80 (3.3% RSD) | match (±4% rounds up from 3.3%) |
| 140 °C, pH 7.1 | "140 °C and 80 psi, holding for 6 s"; pH: control 7.13, processed 7.16 | match |
| "DHS-GC-MS-MS quantitative table" | DHS-GC/MS/MS, MRM, 5-point standard addition (r² 0.96–0.99) | match |

**Bibliography correction:** the canonical citation is Trikusuma, Paravisini & Peterson,
*Food Chemistry* **312**, 126082, **2020** (online Dec 2019, PMID 31901830). The repo's
`2019` suffix is defensible from the online date but the published year is 2020.

**One thing a curator should check:** the values are the **Processed** column (immediately
post-UHT), not the raw control (331 / 59 / 8) and not the 7-week aged sample (683 / 197 /
22). The benchmark's `process_metadata.state: heated_matrix` is consistent with that, so
this looks right — but the aged column differs materially for hexanal and 2-pentylfuran,
so the label is load-bearing and worth an explicit note in the file.

**Residual caveat:** confirmation is against the 2018 thesis, not the 2020 paper (Unpaywall
and OpenAlex both report the paper closed with zero OA locations; ScienceDirect returns 403).
Same authors, same 21-compound DHS-GC/MS/MS design, same UHT conditions, three-for-three on
values *and* on the SD-derived error percentages — the benchmark descends from this data
essentially beyond doubt. Exact fidelity to the published table would need institutional access.

*Access trail:* Crossref (exact record) → Europe PMC `MED/31901830` (abstract only,
"Subscription required") → ScienceDirect 403 → **OhioLINK ETD (success)**. Five OA citing
papers were pulled in full text; none reproduces the concentrations. Semantic Scholar,
ResearchGate (403), OpenAlex and Unpaywall all report no OA copy.

## 2b. Hofmann & Schieberle 1998 — **STILL OPEN**, and the provenance is worse than "undocumented"

`data/benchmarks/cys_ribose_140C_Hofmann1998.json` claims MFT 342 ppb / FFT 200 ppb. The
paper reports **mol %**. The conversion could not be resolved: the body is closed in every
channel checked, and no secondary source in the open literature reproduces its
cysteine/pentose yield table.

**Correct bibliographic record** (Crossref, `10.1021/jf9705983`): Hofmann, T. & Schieberle,
P., "*Quantitative Model Studies on the Effectiveness of Different Precursor Systems in the
Formation of the Intense Food Odorants 2-Furfurylthiol and 2-Methyl-3-furanthiol*",
*J. Agric. Food Chem.* **46**(1), 235–241 (1998).

**The only quantitative text obtainable** — the abstract, via OpenAlex `abstract_inverted_index`:

> "The yields of the two intense food odorants 2-furfurylthiol (FFT) and 2-methyl-3-furanthiol
> (MFT) obtained by heating mixtures of possible precursors in model systems varying in
> temperature, pH value, or water content were determined by using stable isotope dilution
> assays. Although pentoses generated much higher amounts of FFT and MFT than hexoses when
> heated in the presence of cysteine, glucose and rhamnose also gave significant yields.
> Studies on several intermediates indicated the highest yields for MFT (1.4 mol %) when
> hydroxyacetaldehyde and mercapto-2-propanone were reacted for 6 min at 180 degrees C in
> the absence of water. Both intermediates also generated significant amounts of FFT
> (0.05 mol %). However, the system furan-2-aldehyde/H(2)S showed a 10 times higher
> efficiency in generating FFT. Thiamin and norfuraneol/cysteine were less effective
> precursors of MFT."

**The abstract contains no cysteine/ribose number at all.** Its two mol % figures belong to
a *dry* hydroxyacetaldehyde + mercapto-2-propanone system at 180 °C / 6 min — a different
system from the benchmark's.

**The arithmetic, written out.** MFT and FFT are both C₅H₆OS, MW 114.17 g/mol. On the
benchmark's 10 mM (0.010 mol/L) basis, 1 kg/L:

*Reverse — what the repo's ppb imply:*
* MFT: 342 µg/L = 3.42e-4 g/L ÷ 114.17 = 2.9955e-6 mol/L; ÷ 0.010 × 100 = **0.0300 mol %**
* FFT: 200 µg/L = 2.00e-4 g/L ÷ 114.17 = 1.7518e-6 mol/L; ÷ 0.010 × 100 = **0.0175 mol %**

*Forward — what the abstract's numbers would give on the same basis:*
* 1.4 mol % → 0.010 × 0.014 × 114.17 = 1.598e-2 g/L = **15 984 ppb** (46.7x the repo's 342)
* 0.05 mol % → 0.010 × 0.0005 × 114.17 = 5.709e-4 g/L = **571 ppb** (2.85x the repo's 200)

**Verdict: STILL OPEN.** Not confirmed, not refuted. Three things a curator should weigh:

1. **The implied yields are plausible but unsupported.** 0.030 / 0.018 mol % sit ~47x and
   ~29x below the paper's *maximum* (dry, 180 °C, optimised-intermediate) yields. That
   ordering is chemically sensible for an aqueous system — but "not absurd" is all that can
   be said. Nothing traces the numbers to a table.
2. **The condition set appears misattributed.** Cerny 2015 states verbatim: *"in the reaction
   of L-cysteine with molar excess of ribose (molar ratio 1:3), either at 95 °C for 4 h
   (Cerny and Davidek, 2003) or at **145 °C for 20 min (Hofmann and Schieberle, 1995)**, only
   3-mercapto-2-pentanone is produced… In contrast, in the reaction of **equimolar amounts of
   L-cysteine and ribose (at 140 °C for 30 min)**, both isomers are formed **(Mottram and
   Nobrega, 2002)**."* The benchmark's 140 °C / 30 min / equimolar 10 mM is the **Mottram &
   Nobrega 2002** protocol, not Hofmann & Schieberle's 145 °C / 20 min at 1:3. Mottram &
   Nobrega reported headspace GC-MS, not solution ppb, so they are not the source of 342/200 either.
3. **The 10 mM basis is itself unverified.** Since the paper reports mol %, the ppb
   conversion depends entirely on precursor charge and reaction volume, neither of which was
   retrieved. That is a free multiplicative parameter sitting underneath the panel's
   *tightest* contract (1.45x / 0.09 dex).

**Only two clean resolutions:** institutional/ILL access to the body of `10.1021/jf9705983`,
or re-express the benchmark target in the paper's native **mol %** and drop the ppb entirely.
(Reported only — `data/` is owned by the concurrent agent this wave.)

*Access trail:* ACS 403, no SI reachable · Unpaywall `is_oa: false, oa_locations: []` ·
Semantic Scholar `CLOSED` · OpenAlex no repository fulltext · CORE, scholar.archive.org,
colab.ws, x-mol all dead ends · Europe PMC `REF:"jf9705983"` → 5 OA citing papers, all
fetched in full text, **none reproduces the yields** · broader Europe PMC queries
(`"2-furfurylthiol" AND "mol %"`, 26 hits) all false positives · Google Books API, six query
variants, zero volumes · the 1995 and 1997 companion papers are AEDA/FD-factor only ·
`cerny_chapter.txt` cites the paper five times and quotes **no numeric yield** · repo-internal
confirmation that the derivation never existed: `data/lit/slr_incorporation_matrix.json`
entry `hofmann_schieberle_1998_free_mft_fft` records exactly one anchor from this DOI,
`"max MFT 1.4 mol% at 180 C"`, and no ppb value was ever extracted.

## 2c. Brands & van Boekel 2001 — **RECOVERED**. The "unreadable scan" finding was wrong.

The repo's `data/lit/timeseries/brands_sugar_casein_120C_pH68.yml` carries a standing note:

> "The concentration-time multiresponse data that make this thesis valuable … exist ONLY as
> figures in a 1-bit 200 dpi scan. Reading points off those figures would be eyeball
> estimation at low resolution and was NOT attempted."

**That is not true, and the file should be corrected.** The scan *is* bitonal CCITT
(`pdfimages -list`, PDF p.22: `1310 x 1892, gray, 1 bpc, ccitt, 200 ppi`), but rendered at
400 dpi with `pdftoppm` the plot markers are 20–25 px across and unambiguous. I verified
this myself by rendering the page and reading the annotated overlay: the axis, the tick
stubs, the open triangles (glucose) and open squares (fructose) are all cleanly separable,
and the digitizer's marks sit on the marker centres.

**Source correction first.** DOI `10.1021/jf001430b` is *"Reactions of Monosaccharides
during Heating of Sugar−Casein Systems: Building of a Reaction Network Model"* — **not**
"A kinetic modeling study of the Maillard reaction between glucose and casein". **Thesis
Chapter 2 (thesis pp. 9–26 = PDF pp. 16–33) IS that paper, reproduced verbatim**, so the
2001 data were sitting in a file already in the scratchpad. Chapter 4 is the 2002 companion
(`10.1021/jf011164h`, JAFC 50:6725–6739). `brands2002.pdf` in the scratchpad is a
byte-identical duplicate of the thesis, not the 2002 paper.

**Conditions (verbatim, thesis §2.2.2, p.11):** *"Sodium caseinate (3% w/w) and sugar
(150 mM monosaccharide) were dissolved in a phosphate buffer (0.1 M; pH 6.8) to give a molar
ratio of sugar to lysine residues of 10:1. The samples were heated for various times
(0 - 40 min) at 120 °C in an oil bath in screw-capped glass tubes… The reported heating
times include the heating up period of about 2 - 3 minutes."*

### Glucose–casein, 120 °C, initial pH 6.8 (Fig. 2.1), mmol/L

| t (min) | glucose | fructose | formic | acetic | lysine res. | Amadori |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 150.0 | 0.0 | 0.0 | 0.0 | 14.8 | 0.0 |
| 5 | 139.8 | 5.0 | ~0.37 | ~0.09 | 13.67 | 0.67 |
| 10 | 129.9 | 14.2 | ~0.7 | ~0.7 | 12.57 | 1.09 |
| 15 | 122.7 | 19.9 | 0.94 | 1.54 | 11.01 | 1.13 |
| 20 | 110.7 | 23.8 | 1.08 | 2.07 | 10.47 | 0.97 |
| 30 | 103.2 | 29.8 | 2.80 | 3.48 | 8.68 | ~0.79 |
| 40 | 93.8 | 34.0 | 3.25 | 5.23 | 8.53 | ~0.68 |

The fructose–casein table (Fig. 2.2), the pH and total-acid trajectories (Fig. 2.3 —
pH falls 6.73 → 6.45 for glucose and 6.71 → 6.30 for fructose), the A420 browning series
split into total / protein / sugar fractions plus protein-bound melanoidins in mmol/L
(Fig. 2.4), and the isolated-Amadori control (Fig. 2.7) are all in
`<scratchpad>/brands_out_tables.md` (362 lines), with annotated overlay PNGs
(`brands_out_f2*_bbox.png`) so every point can be re-checked by eye.

**Three independent validations that the digitization is sound** (none of them fitted):
1. The t = 0 glucose point reads **150.0 mmol/L** against a stated nominal charge of 150 mM.
2. Fig. 2.3's "total acids by HPLC", digitized off a *different figure with a different
   right-hand axis*, matches formic + acetic from Figs. 2.1/2.2 to **within 0.12 mmol/L at
   all ten time points**.
3. Text cross-check: the thesis states *"a pH decrease of 0.3 pH-unit in a glucose-casein
   system … after heating for 40 minutes"* → 6.73 − 0.30 = 6.43 vs **6.45** digitized;
   *"0.4 pH-unit in a fructose-casein system"* → 6.71 − 0.40 = 6.31 vs **6.30** digitized.
4. Mass balance at glucose 40 min: 93.8 + 34.0 + 8.5 + 0.7 + 5.25 = 142.3 / 150 = **95%**;
   adding the unidentified acids (titration 17.1 − HPLC 8.5) gives ≈101%, matching the
   thesis's own statement of an "almost negligible amount of missing compounds after 40 min".

### The higher-value asset: the fitted kinetic model

Chapter 4's tables OCR'd cleanly and Tables 4.1 and 4.5 plus all eleven Scheme-4.4 ODEs were
verified **digit-for-digit** against 400-dpi page renders. Recovered in full:

* **Scheme 4.1** — 9 ODEs, and **Table 4.1** — 11 rate constants with 95% CIs at 120 °C /
  pH 6.8 (e.g. k1 = 0.01028 ± 0.00046 min⁻¹, k7 = 0.00018 ± 0.00001 L·mmol⁻¹·min⁻¹,
  k8 = 0.11134 ± 0.01146 min⁻¹).
* **Table 4.2** — the full Arrhenius set at **90 / 100 / 110 / 120 / 130 °C** with Ea per
  step (126 ± 2, 117 ± 2, 137 ± 7, 159 ± 5, 129 ± 17, 109 ± 19, 114 ± 2, 124 ± 4, 128 ± 4,
  138 ± 4, 75 ± 11 kJ/mol).
* **Scheme 4.3 / Table 4.4** — the pH-corrected model (pH enters by multiplying every lysine
  term by 10^(−ΔpH)), and **Scheme 4.4 / Table 4.5** — the best model, adding arginine.

That is strictly better than a table of endpoints: with those ODEs and constants the
trajectories can be regenerated on a continuum at any of five temperatures, which makes this
a usable *kinetic* hold-out rather than a handful of points.

**Two defects found in the source, flagged as-is:** Table 4.3 (pH 5.9) is cited in the text
on thesis p.55 but **no such table is printed** on PDF pp. 62–64; and the last two terms of
the printed `d[M]/dt` use `k14` where the rest of the system uses `k15`, almost certainly a
thesis typo. Both are transcribed as printed, not silently repaired.

**Recommended action** (reported only — `data/` is owned by the concurrent agent this wave):
1. Retract the `what_was_deliberately_NOT_extracted` and `extraction.method_note` claims in
   `brands_sugar_casein_120C_pH68.yml`; they assert an impossibility that is false, and they
   are the reason this data sat unused for a wave.
2. Correct `what_a_human_with_access_should_get_next` — it recommends seeking ACS access to
   `10.1021/jf001430b` for vector figures, which is now unnecessary.
3. Promote the trajectories from `brands_out_tables.md` into the timeseries file as
   `data_quality: digitized_figure` (distinct from the existing `table`), carrying the stated
   uncertainties (±0.5 mmol/L on the 0–160 panels, ±0.05 on the 0–8 and 0–16 panels, ±0.01 pH,
   ±0.03 A420; points marked `~` are overlapping markers at ±30% or worse).
4. Consider the Chapter 4 rate constants as a separate artifact. They are a published,
   fully-specified competing kinetic model for exactly the chemistry this repo models, at
   five temperatures — which makes them a far stronger comparator than any single benchmark row.

**Still undigitized, same pipeline applies:** Fig. 4.3 (100 °C out to 120 min), Fig. 4.5
(pH 5.9), Fig. 4.6 (HMF), Fig. 4.7 (pH-stat, 0–240 min), Fig. 4.13 (galactose/tagatose),
Figs. 5.11–5.13 (lactose/maltose).

---

## Mission 2 scorecard

| Target | Outcome |
|---|---|
| **2a Trikusuma 2020** | **CONFIRMED** — 782 / 163 / 24 ppb all match the 2018 OSU thesis Table 6 to rounding; citation year should be 2020, not 2019 |
| **2b Hofmann & Schieberle 1998** | **STILL OPEN** — body closed everywhere; implied yields are MFT 0.0300 mol % / FFT 0.0175 mol % on the assumed 10 mM basis; conditions appear misattributed (they are Mottram & Nobrega 2002's protocol, not Hofmann's 145 °C / 20 min at 1:3) |
| **2c Brands & van Boekel 2001** | **RECOVERED** — full trajectories for both systems plus the complete fitted kinetic model at 5 temperatures; the repo's "unreadable scan" note is false and should be retracted |

---

# ADDENDUM — Wave S1b, 2026-08-27: the routing repair, re-scored

**Everything above this line is the Wave S2 measurement and is left standing verbatim. The
pre-fix numbers are the evidence; deleting them would delete the delta.**

Wave S2's three code findings were all confirmed by inspection and all fixed. No correction
curve was reshaped, no constant refitted, no barrier touched, and the panel was **not**
iterated against — the wiring was chosen from the chemistry, measured once, and reported
whichever way it came out. Two claims flipped to agree; one flipped away; the benchmark
panel got measurably **worse** on four rows and that is reported below rather than
absorbed.

## A1. A correction to the baseline this addendum is measured against

The 18/29 in §2 was scored on `HEAD = c1a12d2`. Wave S1's additive-flux propagator landed
after that, at `HEAD = 263bae8`. **Re-running the identical panel on 263bae8, before any
Wave S1b edit, gives 19/29, not 18/29** — Wave S1 moved SUG-04 (MFT vs FFT in
cysteine/ribose, a *fit-adjacent* row) from agree to `flat`, and moved one independent row
into agreement. So the fit-adjacent bucket had already fallen from 9/9 to 8/9 before this
wave began. **19/29 is the honest baseline for the delta below.** The pH and aw buckets were
untouched by Wave S1 and stood at exactly the 2/7 and 0/3 §3 reports.

## A2. The three defects, each confirmed before it was touched

| # | Wave S2 finding | Confirmed how | Fix |
|---|---|---|---|
| 1 | `get_ph_multiplier` never called on the prediction path | `grep -rn` finds callers only in `kinetics.py`, `pathway_ranker.py`, `cantera_export.py`; none is reachable from `evaluate_benchmark_payload`, which enters `conditions.get_rate_constant()` at `benchmark_validation.py:662` | Called once, from `get_rate_constant`, through a gate (`_enolisation_route_ph_correction`) admitting only the enolisation branch point |
| 2 | `_ionization_correction`'s pyrazine branch unreachable | Enumerated every family the engine emits over all of `data/benchmarks/`: **29 distinct families, not one containing "pyrazine"** | Explicit family set `{Aminoketone_Condensation}` — the family that actually makes 2,5-dimethylpyrazine — at the branch's own unchanged pKa 6.5 |
| 3 | `_water_activity_correction` reaches almost nothing, misses the furan track | Reached 3 of the 29 (`Amadori_Rearrangement`, `Strecker_Degradation`, `Lipid_Strecker_Synergy`); its dehydration branch keyed on `"furfural"`, which matches **no** emitted family either | Membership by **measured net water stoichiometry**: water-releasing families get the (unchanged) `1.3 − aw` dehydration curve, water-consuming families get mass action in `aw`, net-zero families get nothing |

### Which pH term applies where, and why both cannot apply to one family

The two functions are different physics, which is why naive wiring would have double-counted:

* `_ionization_correction` is **reagent availability** — the Henderson–Hasselbalch fraction of
  the nucleophile present as free base. Monotone increasing in pH. Belongs to the families
  where a nitrogen lone pair attacks.
* `get_ph_multiplier` is **route selection** — which way the Amadori compound enolises.
  1,2-enolisation is acid-catalysed and opens 3-deoxyosone → furfural/HMF; 2,3-enolisation is
  base-catalysed and opens 1-deoxyosone → reductone → pyrazine/furanone. A branch-point
  partition, not a reagent count.

The two family sets are disjoint **by construction and by an import-time assertion**, so every
family gets the pH physics exactly once.

`get_ph_multiplier` is deliberately **not** wired to the families its own substrings would
sweep up (`"thiol"`, `"thio"`, `"cysteine"`, `"furan"`, `"condensation"`). The furan track's
acid preference is a property of the branch point; re-applying a 4.9x acid boost at every
downstream thiol addition, cyclodehydration and oxidation would compound one physical effect
five or six times along a single route. And `"condensation"` matches
`Aminoketone_Condensation` — it would have handed the **pyrazine** step the acid-peaked
Schiff Gaussian, i.e. defect 2 again in the opposite direction.

## A3. All dead substring keys found, reported rather than silently deleted

Measured against the 41 family names this engine can emit:

| Key | In | Status |
|---|---|---|
| `"pyrazine"` | `_ionization_correction`, `_water_activity_correction`, `get_ph_multiplier` | **DEAD** — fixed in the first two; documented in the third |
| `"furfural"` | `_water_activity_correction` | **DEAD** — this was the whole dehydration branch |
| `"heyns"` | `get_ph_multiplier` | **DEAD** — no Heyns family is emitted |
| `"nitrogen_heterocycle"` | `get_ph_multiplier` | **DEAD** — never used as a family name |
| `"oxygen_heterocycle"` | `get_ph_multiplier` | **DEAD** — same |
| `"1,2"` and `"2,3"` | `get_ph_multiplier` | **DEAD** — families spell these `1_2` / `2_3`; the underscore forms alongside them do the work, so the comma forms were always decorative |

And one key that is live and was **dangerous**: `"condensation"` in `get_ph_multiplier`'s
Schiff branch matches `Aminoketone_Condensation`. That is the reason the new call site is
gated rather than open.

### A bug this wave introduced and then found, stated because it changes how the numbers read

The first wiring left the Wave H substring guard `_releases_rather_than_attacks_with_the_amine`
in front of `_water_activity_correction`'s new sets. That guard matches `"enolisation"`, so it
returned 1.0 for `Enolisation_1_2` — the 3-deoxyosone → furfural/HMF dehydration, **the furan
track**, the single family finding 3 most needed to reach — before any set was consulted. It
was caught by inspecting the factor table at the snapshot conditions, not by the score. The
guard was removed from that function (the explicit sets subsume it strictly better:
`Enolisation_2_3_Amadori` is net-zero in water, appears in no set, still returns 1.0, and
`tests/unit/test_wave_h_2026_08.py` still pins it). **Every number below is post-fix.**

A second, smaller one: `Furanone_Reductive_Opening` (net +1 water) was omitted from the
water-releasing set by transcription, against the wave's own stated criterion. It was caught
by a test written to re-derive set membership from the enumerated steps rather than trust the
literal in the source, and that test now ships as
`tests/scientific/test_wave_s1b_ph_aw_routing_2026_08.py::test_water_activity_membership_matches_measured_stoichiometry`,
which checks **all 29 emitted families** and currently reports zero mismatches. Correcting the
omission changed **no panel outcome and no benchmark row**; it changed DMHF's aw response
(quoted post-fix below) and nothing else visible.

## A4. THE PANEL, PRE → POST

| Bucket | Wave S2 (`c1a12d2`) | Pre-fix baseline (`263bae8`) | **Post-fix** |
|---|---|---|---|
| **Strictly independent (headline)** | 18/29 (62%) | **19/29 (66%)** | **20/29 (69%)** |
| + system-overlap | 6/6 | 6/6 | 6/6 |
| Screening total | 24/35 | 25/35 | **26/35 (74%)** |
| Fit-adjacent (excluded) | 9/9 | 8/9 | 8/9 |

| Category | Pre-fix | Post-fix | |
|---|---|---|---|
| **ph** | **2/7** | **4/7** | **+2** |
| **moisture_aw** | **0/3** | **0/3** | — |
| **pH + aw combined** | **2/10 (20%)** | **4/10 (40%)** | **+2** |
| sugar_identity | 8/8 | 7/8 | **−1** |
| additive_cysteine | 4/4 | 4/4 | — |
| temperature | 6/8 | 6/8 | — |
| time | 2/2 | 2/2 | — |
| lipid_lane | 2/2 | 2/2 | — |
| matrix_identity | 1/1 | 1/1 | — |
| **Excluding pH and aw** | **17/19 (89%)** | **16/19 (84%)** | **−1** |

**Report these together, as Wave S2 instructed.** The headline rose by one row; the pH bucket
rose by two and the non-pH bucket fell by one. **pH and aw are still worse than chance
(4/10), and aw is still 0/3.** This repair moved the pH failure from *systematically inverted*
to *coin-flip*; it did not make the model a pH advisor. **Wave S2's recommendation to guard
pH and moisture recommendations at runtime stands unchanged.**

### Per-claim flips, every row that moved

| Claim | Observable | Expected | Pre-fix | Post-fix | |
|---|---|---|---|---|---|
| **PH-04** | 2,5-dimethylpyrazine, pH 4.5→5.6→6.5 | increasing | `decreasing` 64.5 / 38.8 / 24.4 | **`increasing`** 0.026 / 0.131 / 2.30 | **GAINED** |
| **PH-06** | 2,5-dimethylpyrazine, pH 4→7→9 | increasing | `decreasing` 99.5 / 16.8 / 14.9 | **`increasing`** 1.80 / 7.95 / 20.7 | **GAINED** |
| **SUG-12** | HMF, fructose vs glucose | A>B | `A>B` 2433 vs 1524 | `A<B` 897 vs 1617 | **LOST** |
| PH-07 | furfural, pH 4→7→9 | flat | `increasing` 752 / 902 / 908 | `peak` 793 / 880 / 792 | still misses, closer |
| AW-02 | acrylamide vs aw | peak | `trough` 8.04 / 3.90 / 4.36 | `decreasing` 3.37 / 0.69 / 0.54 | still misses |
| AW-03 | HMF vs aw | peak | `flat` 570 / 598 / 594 | `increasing` 600 / 636 / 634 | still misses |
| SUG-04 *(fit-adjacent)* | MFT vs FFT, cys/ribose | A>B | `flat` 284 vs 297 | `A<B` 155 vs 267 | already missing, now inverted |

PH-01, PH-02, PH-03, PH-05, AW-01 did not change outcome.

**The two gained rows are exactly the two claims defect 2 was diagnosed against** — the two
independent direct measurements of dimethylpyrazine vs pH (Laemont & Barringer 2023 measure
26.6 → 37.4 → 68.2 ppb over pH 4→7→9). The model now moves them the right way.

### The one lost row, mechanism stated in full

**SUG-12** is HMF from fructose vs glucose at pH 6.0. Glucose reaches HMF through
`Enolisation_1_2`, which is now in the acid-favoured route set and picks up a **3.0x** boost at
pH 6.0. Fructose reaches HMF through `Fructofuranosyl_Dehydration`, a direct ketose
dehydration that bypasses the Amadori branch point entirely and is therefore **excluded** from
that set — so it gets no boost, and glucose overtakes it.

**This is an open owner decision, and the argument runs both ways [P].**
*For inclusion:* fructofuranosyl dehydration to HMF is the textbook acid-catalysed hexose
dehydration; two routes to the *same product* now receive opposite pH treatment, which is an
internal inconsistency. *For exclusion (what shipped):* `get_ph_multiplier`'s branch is
documented as the 1,2- vs 2,3-**enolisation** partition, and this family is not part of that
partition. **The exclusion was written down before the panel was re-scored, and was
deliberately not revisited afterwards** — re-wiring it now, having seen that it costs a row,
would be exactly the optimisation-against-the-panel this document exists to prevent. The
owner should decide it on the chemistry.

### Why aw is still 0/3 — a structural finding, not a routing one

AW-01 and AW-03 both score **HMF** in a glycine/glucose system. The correction now reaches
that system's chemistry properly (5 of its 7 emitted families, up from 2), and it visibly
moves compounds: DMHF falls **115.4 → 45.5 → 49.6 ppb** over aw 0.25 → 0.65 → 0.95 — a
2.5x separation where there was none — and dimethylpyrazine falls 2.15 → 0.29 → 0.26. But
**HMF cannot respond**, because HMF and
furfural are both products of `Enolisation_1_2` — the *same* family — so they always carry the
*same* aw and pH factor, and together they are **90–96%** of that system's volatile budget.
Their shares are pinned against each other, and the projection budget itself is
aw-independent.

**No family-level correction can move HMF vs aw in a two-precursor system.** That is a
statement about the allocation layer, not about the aw physics, and it is why AW-01/AW-03 were
not going to be fixed by this wave regardless of how the families were assigned. Recorded as
[P]: the budget has no moisture dependence.

AW-02 (acrylamide vs aw) is a separate matter — it needs the acrylamide **destruction** term
Wave S2 identified as the model's most consequential single miss, and this wave did not add one.

## A5. THE COST — the benchmark panel got worse, on four rows

Regenerated `benchmark_summary.json`. **`status_counts` are unchanged (scale-gap 8 /
pass-no-ranking 2 / pass 4), strict-ready is unchanged at 0/14, and the four internal
snapshots recover exactly after the documented refresh.** Four real rows moved, all four in
the wrong direction:

| Benchmark | max_ratio | MALE | |
|---|---|---|---|
| `cys_ribose_140C_Hofmann1998` | 1.4864 → **2.2086** | 0.1267 → **0.2352** | WORSE |
| `thiamine_cys_glucose_120C_Bolton1994` | 748.02 → **6730.85** | 2.8739 → **3.8281** | MUCH WORSE |
| `thiamine_cys_xylose_145C_Cerny2008` | 2.7874 → **23.4061** | 0.4452 → **1.3693** | MUCH WORSE |
| `resconi_2023_pbma_beef_identity_benchmark` | 4.4024 → **5.4570** | 0.6437 → **0.7370** | WORSE |

**Nothing was refitted to claw any of this back.** The Hofmann contract (1.45x / 0.09 dex) was
already failing on both criteria after Wave S1 and fails by more now; it was not relaxed.

> **Addendum, 2026-08-27 (Wave S2c) — the Hofmann row's yardstick was retired the same day.**
> Wave S2b established that `cys_ribose_140C_Hofmann1998`'s MFT 342 / FFT 200 ppb are a
> repo-internal derivation: interior points of two invented, overlapping mol % bands in
> `data/benchmarks/maillard_validation_benchmarks.md` §1.3, an abstract-reconstructed table
> committed in the same commit as the benchmark JSON. Wave S2c retired its 1.45x / 0.09 dex
> contract (it was ~1.7x *tighter* than the 2.5x spread of the band its own target came from),
> demoted its tier `PRIMARY → REFERENCE`, marked both values `no_verifiable_source`, and
> reverted `thiol_addition_pentodiulose` 26.35 → 28.60 because that constant's sole fit target
> was this benchmark. On the post-revert tree the Hofmann row reads **4.3797 / 0.4041 dex** and
> the Cerny row **25.741 / 1.4114** — both worse again, both published. **Nothing in Wave S1b's
> attribution changes:** the pH/aw routing repair really did move those rows for the reasons
> given here. What changes is that the Hofmann row cannot be read as the model disagreeing with
> a measurement, in this table or anywhere else. **The sulfur branch has zero absolute
> literature anchors.**

### Attribution, measured one fix at a time

Each fix was measured alone by emptying the other two family sets at runtime (no source edit):

| System / observable | none | fix 1 only | fix 2 only | fix 3 only | all three |
|---|---|---|---|---|---|
| Hofmann MFT (meas 342) | 283.6 | 219.7 | 296.5 | 286.2 | **154.8** |
| Hofmann FFT (meas 200) | 297.3 | 317.6 | 310.9 | 245.2 | **267.5** |
| Bolton MFT (meas 13.0) | 0.01738 | 0.01659 | 0.01744 | **0.002055** | **0.001931** |
| Resconi furfural | 3149 | 3721 | 3265 | 2922 | **3903** |

* **Bolton and Cerny are carried almost entirely by fix 3 (aw).** At aw 0.98 every
  water-shedding step is multiplied by `1.3 − 0.98 = 0.32`, i.e. **+0.89 kcal/mol** of
  effective barrier, and the thiamine → MFT lane's terminal `Furan_Ring_Aromatisation` takes
  it while the `Additive_Thermal_Degradation` steps upstream of it do not (that family is
  stoichiometrically non-uniform — +2/0/−1/−2 across its steps — so it is excluded, since one
  family-level factor cannot honestly represent it). MFT is a small share of a large budget in
  that system, so its competitors absorb what it loses and the ratio amplifies.
* **Hofmann is carried by fix 1 (route pH), with interaction.** At pH 5.0 `Enolisation_1_2`
  gets a **4.5x** acid boost and `Enolisation_2_3_Amadori` gets ≈1.0, so the FFT/furfural arm
  gains on the MFT arm, against a benchmark whose targets are MFT 342 > FFT 200 at pH 5.0.

  **CORRECTION, added after Wave S2b (2026-08-27, same day).** This bullet originally called
  that "the textbook enolisation pH physics disagreeing with the Hofmann benchmark — one of
  the two is wrong". That framed it as a conflict with a measurement, and **there is no
  measurement on the other side.** Wave S2b traced 342 / 200 ppb to this repository's own
  `data/benchmarks/maillard_validation_benchmarks.md` §1.3 — an abstract-reconstructed range
  table (MFT `~0.02–0.05` mol %, FFT `~0.01–0.03` mol %) committed in the *same commit* that
  created the benchmark file. Both shipped values are interior points of two **overlapping**
  invented bands (MFT 228–571 ppb, FFT 114–342 ppb), so the MFT > FFT ordering is an artifact
  of midpoint selection, and the benchmark's own 1.45× / 0.09 dex contract is ~1.7× tighter
  than the band it came from. The mechanism above is real and the degradation is real and
  reported; what is not available is a measurement saying the pH physics is wrong. Wave S2b
  also corrects §2b of this report by one paper: Cerny 2015's "145 °C / 20 min at 1:3"
  sentence cites Hofmann & Schieberle **1995** (`10.1021/jf00056a042`), not the 1998 paper, so
  Cerny says nothing about this benchmark's conditions. See `## Wave S2b` in
  `tasks/audit_remediation.md` and the `content_verification_note.wave_s2_followup` block in
  the benchmark file, which carries both corrections.

### The internal snapshots

They are model-generated reproducibility baselines, refreshed per the Wave M sequence
(`GENERATOR_TAG` v9 → **v10**). Movement, pea and soy identically:

```
2,5-dimethylpyrazine              x0.0022206     <- penalised twice: 89% protonated at pH 5.6 AND sheds 2 H2O
2-methyl-3-furanthiol             x0.337600
bis(2-methyl-3-furyl) disulfide   x0.496260
2-furfurylthiol                   x0.504180
furfural                          x1.074000      <- its 1,2-enolisation route gets the acid boost
Hexanal, Nonanal                  x1 (unchanged) <- injected by the lipid lane, never propagated
```

**The ranking contract changed, and that is a real movement rather than a refresh artifact:**
2-furfurylthiol now outranks 2-methyl-3-furanthiol, and the disulfide now outranks
2,5-dimethylpyrazine.

## A6. What this addendum does and does not change about §7's licence

Unchanged: it licenses **screening** claims of the form "which arm gives more of X" with pH and
moisture held fixed, and **nothing quantitative**. What changed is narrow and should be stated
narrowly: **the pH direction for nitrogen heterocycles is no longer systematically inverted.**
pH is 4/7, aw is 0/3, and 4/10 combined is still at or below chance. **The model still licenses
no pH or moisture recommendation.** The runtime guard Wave S2 asked for is still the right call.

---

# CURRENT STANDING — the panel as it scores on the shipped tree

**Added 2026-08-28 (Wave S5).** Everything above this line is left verbatim: §2 is the Wave S2
measurement, the addendum is the Wave S1b re-score. Neither was ever updated for the wave that
landed after them, so a reader who stopped at the addendum carried away **20/29** when the
shipped tree scores **21/29**. This section is the single place that states the *current* number,
and it is the section machine-read by `src/directional_reliability.py`, which is what the
`maillard` CLI prints as its per-axis reliability tags. **Nothing here is a new measurement.**
The one row that moved is the Wave T4 Heyns-consistency fix, recorded with its mechanism in
`tasks/audit_remediation.md` (Wave T4): `SUG-12` (HMF, fructose vs glucose) flipped MISS → OK,
carrying `sugar_identity` 7/8 → 8/8 and the headline 20/29 → 21/29. Waves S3 and S4 both
re-scored the panel and both found it **byte-identical** at 21/29.

**Updated 2026-08-28 (Wave W).** Eleven claims were added to the panel from five
interlibrary-loan papers whose full text this repository had never held (`HOF-01`…`HOF-04`
from Hofmann & Schieberle 1998, `MOT-01`…`MOT-06` from Mottram & Nobrega 2002 — see §W below
for what each one is and why four of them are unevaluable). Six are evaluable and all six are
strictly independent. **They scored 3/6, which is worse than the panel average and pulls the
headline down from 21/29 (72%) to 24/35 (69%).** The movement is not evenly spread and the
direction is the interesting part: the two new pH rows both AGREE, taking `ph` from 4/7 to
6/9, while `sugar_identity` — the category this report called *Trustworthy* at 8/8 — takes its
first two misses and drops to 9/11. **The category the panel was most confident about is the
one the new evidence damaged, and the category it called broken is the one that improved.**

| bucket | agree | evaluable |
|---|---:|---:|
| strictly independent (headline) | 24 | 35 |
| system-overlap | 6 | 6 |
| fit-adjacent (excluded from the headline) | 8 | 9 |
| independent, excluding ph and moisture_aw | 18 | 23 |
| independent, ph and moisture_aw only | 6 | 12 |

**Read the last two rows, and the category table below, with their denominators in mind — they
are different denominators and conflating them inflates the score.** The
category table sums to **30/41**, which is the *screening* bucket: the 35 independent claims
plus the 6 system-overlap ones. Excluding pH and water activity it sums to 24/29 (83%). But
all six overlap rows sit outside pH and aw, so the strictly-independent version of that same
statement is **18/23 (78%)**, which is the number to quote alongside the 24/35 headline. All
twelve pH and aw rows are independent, which is why 23 + 12 = 35 exactly.

| category | agree | evaluable |
|---|---:|---:|
| sugar_identity | 9 | 11 |
| additive_cysteine | 4 | 4 |
| temperature | 6 | 8 |
| time | 2 | 2 |
| lipid_lane | 2 | 2 |
| matrix_identity | 1 | 1 |
| ranking | 0 | 1 |
| ph | 6 | 9 |
| moisture_aw | 0 | 3 |

Stated in one line, on the independent bucket throughout: **24/35 overall; 18/23 (78%)
excluding pH and water activity; 6/12 (50%) on pH and water activity, which is at chance;
0/3 on water activity alone.** Twelve further claims are unevaluable and are counted in
neither column (§7, §W).

**What moved, and it is not flattering.** Before Wave W this section read *21/29 overall;
17/19 (89%) excluding pH and water activity*. The 89% was the number this report leaned on
hardest — it is the basis of the §8 licence to use the model for screening with pH and
moisture held fixed. Six genuinely out-of-sample rows took it to 78%. The `sugar_identity`
category, previously **8/8 and labelled *Trustworthy* without qualification**, is now 9/11:
the model separates glucose from fructose by 52× (MFT) and 154× (FFT) where Hofmann &
Schieberle measure them within 1.3× of each other, in the same table, the same run and the
same lab. That is not a near miss on an ordering; it is two orders of magnitude of spurious
structure between two hexoses. The `ranking` category appears for the first time with a
single independent row and it is a miss (`MOT-03`). Against that, `ph` improved from 4/7 to
6/9 on two rows that AGREE — the first pH rows this panel has ever scored from a source with
a stated buffer system, and they agree with the Hofmann pH hold-out's direction as well.

**The licence in §8 and §A6 is unchanged by this section.** The model licenses screening claims
of the form "which arm gives more of X" with pH and moisture held fixed, and nothing
quantitative. It licenses no pH recommendation and no moisture recommendation.

### How to update this table

Re-score the panel with the scratchpad runner described in §9, then edit the two tables above
and say in one line which rows moved and why. The verdict thresholds the CLI applies to these
counts (`trust` / `caution` / `do-not-use`) live in `src/directional_reliability.py` and are
stated in the README model card; they are deliberately *not* stored here, so that this file
stays a record of measurement and nothing else.

---

## §W. The Wave W rows (2026-08-28)

Eleven claims added from five papers the repository obtained by interlibrary loan on
2026-08-28 and read in full for the first time. **Nothing in the model was tuned to build
these rows, and Wave W fitted no constant to anything.** All eleven are `fit_status:
independent`.

| Claim | Category | Observable | Literature | Model | Outcome | Predicted values (ppb) |
|---|---|---|---|---|---|---|
| HOF-01 | sugar_identity | MFT | `A>B` | `A<B` | **MISS** | cysteine + D-fructose = 1.729; cysteine + D-glucose = 90.23 |
| HOF-02 | sugar_identity | FFT | `A>B` | `A<B` | **MISS** | cysteine + D-fructose = 5.388; cysteine + D-glucose = 828.3 |
| HOF-03 | sugar_identity | FFT | `A>B` | `A>B` | agree | ribose dry 180 °C = 5545; glucose dry 180 °C = 2897 |
| HOF-04 | moisture_aw | MFT, FFT | `A>B` | `-` | n/a | — (triply confounded in the source; see below) |
| MOT-01 | ph | MFT | `A>B` | `A>B` | agree | pH 4.2 = 1258; pH 5.6 = 103.2 |
| MOT-02 | ph | FFT | `A>B` | `A>B` | agree | pH 4.2 = 1720; pH 5.6 = 454.8 |
| MOT-03 | ranking | MFT | `A>B` | `A<B` | **MISS** | MFT = 103.2; FFT = 454.8 |
| MOT-04 | sugar_identity | MFT | `A>B` | `-` | n/a | — (ribose-5-phosphate not in the resolver) |
| MOT-05 | sugar_identity | MFT | `A>B` | `-` | n/a | — (IMP not in the resolver) |
| MOT-06 | ph | MFT | `A<B` | `-` | n/a | — (no buffer term; no pH trajectory) |

**3 agree / 6 evaluable. Four unevaluable, and each one names what would make it evaluable.**

### The two misses that matter most

**HOF-01 / HOF-02 — the hexose pair.** Hofmann & Schieberle's Table 1 measures fructose and
glucose with cysteine under identical conditions: FFT 3.2 vs 2.8 µg, MFT 2.5 vs 1.9 µg. Both
pairs are within 1.32×, which is at the edge of the paper's own ±10% replicate spread — the
honest reading of the source is *these two hexoses are the same to well within a factor of
two*, with fructose marginally ahead. The model puts glucose ahead of fructose by **52× on
MFT and 154× on FFT**. The sign is wrong and the magnitude is wrong by two orders. This is
the first evidence the panel has ever held on hexose-vs-hexose discrimination, and the model
fails it badly. It also demotes the previous `sugar_identity` 8/8: every one of those eight
rows was pentose-vs-hexose or thiol-vs-thiol, i.e. comparisons across a large real gap. The
model is good at large gaps and inventing large gaps where none exist is the same defect
seen from the other side.

**MOT-03 — the MFT/FFT branch ratio.** Both Hofmann (SIDA, absolute) and Mottram (Tenax,
relative) put MFT above FFT in buffered cysteine/ribose near pH 5; the model puts FFT above
MFT by 4.4×. Note the caveat carried on the claim itself: Mottram's response factor of 1
makes his cross-compound comparison weak evidence. Hofmann's SIDA numbers are the authority
here (19.8 vs 12.1 µg at pH 5.0), and they agree with Mottram's sign, so the model is
contradicted by both. The corresponding absolute rows are scored in
`data/benchmarks/hofmann1998_ribose_cysteine_145C_20min_pH5.json`, where the model
over-predicts FFT by 12.3× and MFT by 2.4× — the same defect measured on the absolute scale.

### The unevaluable four are the finding, not a gap in the work

- **HOF-04** (dry vs aqueous within one sugar) changes temperature, time *and* water content
  at once, and the source assigns no water activity to the dry arm. A hit would be luck and
  a miss uninterpretable. *Would become evaluable with:* a dry-heat source that reports an
  a<sub>w</sub>, or an aqueous arm at matched 180 °C / 6 min.
- **MOT-04 / MOT-05** need ribose-5-phosphate and inosine-5′-monophosphate, neither of which
  is among the 19 species `src/precursor_resolver.py` resolves. MOT-05 is the sharper one:
  IMP is the dominant ribose reservoir in post-mortem muscle, so the model cannot reason
  about the actual precursor pool of cooked meat. *Would become evaluable with:* the species
  **and** a dephosphorylation step — the species alone would produce a confidently wrong
  answer, since the paper's own mechanism is dephosphorylation.
- **MOT-06** (buffer catalysis) needs both a buffer term the model does not have and a pH
  *trajectory* it cannot represent: the unbuffered arm drifts from pH 5.6 to 3.9 during the
  run and the model integrates at one fixed pH. No pH was assigned to it, deliberately —
  picking 5.6 or 3.9 would have been this repository inventing the number it then scored
  itself against.

### Why Mottram & Nobrega is ordinal-only, permanently

Wave S2b predicted at ~85% confidence that this paper's values were Tenax-relative. The paper
settles it in its own words: *"This allowed comparison of the relative contributions the
volatiles made to the headspaces of the different systems but did not provide absolute
concentrations in the aqueous solutions."* It can never become an absolute benchmark here, at
any future date — the quantity it reports is not the quantity the model predicts. See
`sources.mottram2002.ordinal_only_reason` in the panel YAML for the full quote.
