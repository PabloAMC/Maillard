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
