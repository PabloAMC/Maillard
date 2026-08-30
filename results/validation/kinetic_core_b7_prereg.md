# Build Wave B7 — THE FURANIC CHANNEL — PRE-REGISTRATION

**Written 2026-08-29, BEFORE any fit was run and before any hold-out was scored.**
**Branch `audit-remediation`, HEAD `49df685`.**

Module 6 of the kinetic-core rebuild: the **HMF node** (K5a) and the **DMHF/furanone node**
(K5b). Roles are those of `docs/reference/FIT_HOLDOUT_DECLARATION.md` **Amendments 8 and 12**;
the evidence is `data/lit/extraction_dossiers/k5a_hmf_synthesis.md` and
`k5b_dmhf_synthesis.md` plus the seven per-paper dossiers they consolidate.

---

## §0. EXPOSURE DISCLOSURE — read this first

**§0.1 — `results/validation/cutover_final_exam.md` was opened while locating the seven
refused bundles, and it prints the MEASURED value of every one of the 40 exam points,
including all 5 HMF rows and both DMHF rows.** The instruction was never to open
`data/benchmarks/external_validation/` outside the exam generator, and that was obeyed —
no bundle file was read. But the shipped exam ARTEFACT carries the same numbers in its
"Every point, old lane vs core" table, and this builder saw them before writing a line of
`src/`.

This is the same class of exposure as **Amendment 10 clause 1** (Frankel's tocopherol
columns), and it is handled the same way: **structurally, not procedurally.**

* **The HMF node has NO fitted parameter at all.** Its seven rate constants are INGESTED
  WHOLE from Kocadagli & Gokmen 2016's glucose Table 2 and Hamzalioglu & Gokmen 2018's
  Table 1 refit. There is no free number to tune toward 12 280 / 57 250 / 101 000 / 7 000 /
  17 400 ppb even in principle.
* **The DMHF node has exactly ONE fitted parameter**, `k_dpo_af`, and it is calibrated on
  three cells of Blank 1997 Table 1 in `ug per mmol of PENTOSE` at 90 °C — a system with no
  glucose, no alanine, no asparagine and no relation to any exam bundle. Its objective is
  written out in §2 below, in full, before it was run.
* **A literal-grep firewall test** (`tests/unit/test_kinetic_core_b7.py`) asserts that no
  exam-side measured literal appears anywhere in the furanic package or in the frozen fit
  report, on the B2.1–B2.4 pattern.

**Under the Amendment 9 clause 1 / Amendment 10 clause 1 precedent, the seven exam rows are
recorded as `seen_diagnostic` and may never gate.** They are still reported, per bundle,
because the wave was asked for that delta — but the conclusion this wave is entitled to draw
from them is "the refusals became answers/sharper refusals", not "the module passes".

**§0.2 — Every hold-out below was ALSO printed in a dossier this builder was instructed to
read.** K5a §9.1 and K5b §8.3 name each hold-out with its value attached. That is unavoidable
for a wave whose brief is "read these syntheses, then build". The mitigations are the same
two: (i) nothing in the module is fitted to anything but the three Blank cells, so there is
one degree of freedom in the entire channel and it is spent on a pentose system that appears
in none of the hold-outs; (ii) each hold-out below carries a **written, quantitative
prediction made before the scorer existed**, so the comparison is against a commitment rather
than an impression.

---

## §1. WHAT IS BEING BUILT

### 1.1 The HMF node — two parallel sources, one sink

```
  Glc ──k_glc_tdg──► 3-DG ──k_tdg_ddg──► 3,4-DG ──k_ddg_hmf──►  HMF
                      │      (RDS, C3)            (fast, Ea=0    ▲
                      │                            AUTHOR-FIXED)  │
                (also from the Amadori compound,                  │
                 B1's r_ama_tdg — TWO parents, no fixed split)    │
                                                                  │
  Fru ──k_fru_int──►  Int  ──────────k_int_hmf────────────────────┘
        (fast, C4)  (UNMEASURED — C2)   (RDS of this limb)

  sink:  HMF + Cys ──k_hmf_cys──► adduct     (2nd order, CLAMPED at 50 °C)
         HMF       ──k_hmf_self─► products   (Ea = 0, single-temperature)
```

Adopted **without modification** from the source topology that **four independent groups in
six matrices all write** (K5a §8.1). **There is no branch-fraction constant and there cannot
be one**: the share each limb takes is whatever the dynamic `Fru` and `TDG` pools make it,
which is what every paper that explains its own verdict actually appeals to.

**Every HMF formation constant is INGESTED, none is fitted.** Kocadagli glucose Table 2
steps 3, 4, 5, 6, 7, 8, 11, re-referenced from T_b = 180 °C to the core's 100 °C once, in
`parameters_furanic._at_100C`.

### 1.2 The DMHF node — three edges

```
  EDGE A (intact skeleton, via acetylformoin)
     hexose:  1-DG ──k_odg_af──► acetylformoin ──k_af_dmhf──► DMHF
     pentose: 1-DPO + Gly ──k_dpo_af──► acetylformoin ──► DMHF     [C5 + C1 → C6]
  EDGE B (methylglyoxal, C3 + C3, acetylformoin-free BY CONSTRUCTION)
     2 MGO ──k_mgo_dmhf──► DMHF
  EDGE C (cysteine / H2S sink)
     DMHF + H2S ──k_dmhf_h2s = 0.0 EXACTLY──► thiophenone
```

* **`k_dpo_af` is the only fitted number in the channel.**
* **`k_odg_af` is DERIVED** as `k_dpo_af × 1000 mmol/L` — the pentose edge's
  pseudo-first-order constant at Blank's 1 M amine loading. A declared transfer, because
  **no absolute hexose DMHF yield exists in any of the five papers.**
* **`k_mgo_dmhf` is a DIGITISED PRIOR** from Wang & Ho's Figure 1 bar chart, carried with a
  flag and never fitted.
* **`k_dmhf_h2s` is EXACTLY ZERO.** Fitting anything to Shu & Ho's 6.0 % GC area is a
  prohibited derivation named by name in K5b §8.6.
* **No furanone edge carries an activation energy from any source.** All five papers are
  single-temperature. The partition Ea is **0 by declaration** and DMHF's whole temperature
  dependence is inherited from the MEASURED 1-deoxyosone-formation steps. The assumption is
  banded at ±50 kJ/mol and **priced by re-integration**, exactly as B6 prices its Q10.

### 1.3 What this wave will perturb, and it is not zero

The furanic block hangs on the **trunk**, so its eleven steps run inside the trunk, the
sulfur and the acrylamide lanes. Four of them put a new drain on a B1-fitted species:
`Fru → Int` (~12 % of fructose's flux at 100 °C, ~15 % at 180 °C), `3-DG → 3,4-DG` (~11 % of
3-DG's sink), `3-DG → MGO` (~7 %) and `Glc → 3-DG` (~3 %). **B1 is NOT refit** — its four
constants are frozen and this wave has no licence to move them — so these appear as small
movements in every previously-scored trunk/sulfur/acrylamide number. **Pre-registered
expectation: no previously-answered exam row moves by more than 1.5× on this account.** The
exam is reported both ways where it does.

---

## §2. THE FIT — declared BEFORE it ran

### 2.1 The objective

**FIT ROWS (three, and only three):** Blank, Fay, Lakner & Schlosser 1997, JAFC 45:2642–2648,
Table 1, the **pentose + glycine HDMF** cells — arabinose **5.1**, xylose **2.6**, ribose
**3.6** µg per mmol of sugar. Isotope-dilution assay against synthesised [¹³C₂]HDMF,
stated maximum SD ≤ 10 %, n ≥ 2 assays × 2 injections. 90 °C, 1 h, 0.2 M phosphate, pH 6,
1 M sugar + 1 M glycine.

**OBJECTIVE:** sum of squared log₁₀ residuals, `σ_log10 = log10(1.10) = 0.04139` from the
paper's own stated maximum SD. One free parameter, `k_dpo_af`. **Two seeded starts**
(`1e-8` and `1e-4` L mmol⁻¹ min⁻¹, seed 20260829), Nelder–Mead on `log10 k`, and the two
starts must agree to 1e-6 in `log10 k` or the fit is reported as non-identified.

**B2.4's weighting convention** (`PH_ENDPOINT_WEIGHT` and the pH-unit/log-fold exchange rate)
is **not engaged**: every fit row here is in one unit (µg/mmol) and none is a pH endpoint, so
there is no mixed-unit exchange to declare. Stated explicitly so its absence is a decision
rather than an omission.

### 2.2 What is deliberately NOT fitted, though it is a declared-FIT row

Nine of Blank 1997 Table 1's twelve cells are **unrepresentable** by the core and are
reported as such rather than dropped:

| cells | why the core cannot represent them |
|---|---|
| all six **alanine** cells | pentose lives on the sulfur lane and alanine only on the acrylamide lane; the two do not compose. A lane-algebra limitation, not a chemistry one. |
| the three **glycine EHMF** cells | HEMF needs a C2 Strecker donor, i.e. alanine — same lane conflict. |
| the **1.96× sugar spread** inside the three fitted cells | the sulfur lane carries ONE generic aldopentose. Reported as a fit residual, exactly as B2 reports Hofmann's 1.38× ribose/xylose gap. |

**Pre-registered consequence:** the fit cannot do better than the geometric mean of 5.1 / 2.6
/ 3.6, so its three residuals will be **≈ +0.16, −0.13 and 0.00 decades** and its RMS will be
**≈ 0.12 decades ≈ 1.33×**, entirely from the sugar-identity axis the model does not carry.
**If the RMS comes out materially below that, something has been fitted that should not have
been.**

### 2.3 The unconstrained constant, and its sweep

`k_af_dmhf` (acetylformoin → DMHF) has **no source of any kind**. It encodes "acetylformoin
does not accumulate" and is set to 10× the parent deoxyosone's own measured sink.
**Pre-registered check:** sweeping it over **three decades** (0.1× to 100× that value) must
move the DMHF prediction by **less than 1.2×**. If it moves more, the constant is
rate-limiting after all, the assumption is false, and the module must say so.

---

## §3. THE HOLD-OUTS — predictions written before the scorer existed

All seven are **blind**: no scorer for any of them existed when this section was written.

### H1 — Kocadagli's glucose-**NaCl** arm (the single-variable perturbation)

**The model has no ionic-strength term, no salt species and no activity coefficient.**
It therefore predicts, by construction and with no free parameter:

* `k(Fru→Int)` NaCl/glucose ratio = **exactly 1.00** at 160, 180 and 200 °C;
* HMF mole-conversion NaCl/glucose ratio = **exactly 1.00** at all three temperatures.

**PRE-REGISTERED OUTCOME: 0 of 3 rate-ratio cells inside the 3× band, and 2 of 3
mole-conversion cells inside it — the 180 °C and 200 °C ones, and by arithmetic rather than
by chemistry** (the measured NaCl enhancement of the HMF yield collapses as temperature
rises, so at the top of the ladder "no effect" is nearly right for the wrong reason). This is
a **declared, quantified FAIL**, pre-registered as such, and it is exactly what a hold-out on
a variable the model does not carry is for.

### H2 — Gursul Aktag 2020's **27 °C zero-accumulation** row

The cheapest and most informative single test in the module: a model with an over-large HMF
activation energy fails it immediately. Fruit-juice composition, pH 3.4, **27 °C, 24 weeks**
— 120–170 °C below every fit-panel temperature, and the module's only genuine temperature
extrapolation.

**PRE-REGISTERED PREDICTION: PASS.** Arrhenius arithmetic on the ingested constants
(Ea 100.4 on `Fru→Int`, 151.4 on `Int→HMF`) gives HMF **below 100 ppb** at 27 °C / 24 wk,
against a measured *no detectable accumulation*. The declared pass floor is **100 µg/L**,
set here and not adjusted later. ⚠️ Gursul Aktag's own LOD/LOQ pair is **internally
impossible** (an LOD three orders above its own LOQ), so the floor is this module's and is
declared rather than taken from the paper.

### H3 — Hamzalioglu's **matrix-vs-water selectivity pair** (11.4× → 1.2×)

A same-method matrix/water pair on a rate ratio — the class `k3` §C.2 says the corpus lacks.

**PRE-REGISTERED PREDICTION: FAIL, by construction, with the missing ingredient named.**
The model carries **one** amine on the sulfur lane (cysteine), so a Cys : Arg : Lys
selectivity cannot be formed at all; and it carries **no moisture, matrix or water-activity
term** on the sink edge, so `k_Cys(coffee)/k_Cys(water)` is predicted at **exactly 1.000**
against a measured collapse. The scorer reports the predicted 1.000 beside the measured ratio
and the fold error between them.

### H4 — the **shu1988 × wang2008 paired sink/net-pH** test

Two experiments from opposite sides of the node, two decades apart, same lab: Shu's sink
markers go OFF by pH 7.1, and Wang & Ho's **net** DMHF in a cysteine system RISES with pH.
Scored as **ONE paired test**; a model that reproduces one and not the other has the sink or
the formation edge mis-signed, and the pair says which.

**PRE-REGISTERED PREDICTION:**
* **sink half — REFUSED.** Edge C's rate is exactly zero, so the model makes no prediction
  about the pH shape of a sink it does not run. This is a *sharper refusal* than B6's: the
  edge, the species and the balanced stoichiometry all exist, and what is missing is named
  (a magnitude, from Haleva-Toledo 1999).
* **formation half — the model will predict net DMHF RISING with pH, i.e. the RIGHT SIGN FOR
  THE WRONG REASON**, and this wave commits to saying so. The rise will come from the sulfur
  lane's **base-favoured 2,3-enolisation** pH factor on `r_pent_dpo`/`r_arp_dpo`, which is
  Edge A. Wang & Ho's curve is Edge B, in a system with no sugar at all. **Sign agreement is
  therefore recorded as coincidental and NOT as a pass.**

### H5 — `apriyantono1993_xylose_lysine_pH_trajectory_pair` (the named role)

The corpus's **only** held-vs-drifting pair, scored as **ONE paired log-ratio test** across
four channels. It also exams the B2.2/B2.3 **pH state**, because it asks whether the model's
internal pH EVOLVES correctly rather than whether it gets two buffered points right.

**Declared limitation, stated before scoring:** the paper's amine is **lysine**, which lives
only on the acrylamide lane while the pentose lives only on the sulfur lane. **The pair is
run sugar-only**, and the amine's contribution is an uncontrolled, declared difference. Only
direction and order of magnitude are scored.

**PRE-REGISTERED PREDICTIONS, per channel:**

| channel | measured | model | prediction |
|---|---|---|---|
| total volatiles | 143× UP on drift | **REFUSED** — no total-volatile observable exists | no score |
| **2-furfural** | 274× UP on drift | `r_pent_tdp` carries the **acid** pH factor | **DIRECTION PASS, MAGNITUDE UNDER-SHOOT.** Predicted ratio > 1 but ≪ 274× |
| N-heterocycles (16 cpds) | 75× DOWN | **REFUSED** — no pyrazine, pyrrole, pyridine or pyrrolizine species | no score |
| **norfuraneol** | present → not detected | `r_pent_dpo` carries the **base** factor | **DIRECTION PASS** (down on drift), magnitude unscorable (measured value is at the floor) |
| **the pH state itself** | held 5.0 vs drifting 4.9 → **2.6** | B2.2 dynamic pH, unbuffered arm | **DIRECTION PASS, MAGNITUDE UNDER-SHOOT** — the model will drift far less than 2.3 pH units |

**Scored as 2 of 4 chemical channels plus the pH-trajectory channel; 2 refused with the
missing species named.**

### H6 — **norfuraneol ≫ DMHF** at the deoxypentosone fork

Ordinal only, forever: supported by two papers and six systems and **quantified in neither**.

**PRE-REGISTERED PREDICTION: PASS, by a large margin.** The calibration puts DMHF at
~10⁻⁵ mol % of the pentose charge while norfuraneol is a major product of the same
deoxypentosone. ⚠️ **Apriyantono's norfuraneol cell must NOT be scored against this** (K5b
§8.5) — both terms are at the detection floor, the amine is lysine, and the isolation is SDE.

### H7 — the **cutover exam re-run**: 5 HMF + 2 DMHF bundles

| bundle | compound | pre-registered outcome |
|---|---|---|
| `mp_holdout_glucose_only_autoclave_121C_Steinha` | HMF | **ANSWERED.** The amine-free case the FIT panel is literally made of; the closest match in the whole exam. Expect the smallest fold error of the five. |
| `mp_holdout_glucose_asparagine_180C_30min_water` | HMF | **ANSWERED.** Asparagine drives the acrylamide lane; HMF comes from the amine-free limbs. |
| `mp_holdout_fructose_asparagine_180C_Lin2022` | HMF | **ANSWERED**, and expected **HIGHER** than the glucose systems: fructose feeds the dominant limb directly, and Agcam measures a 4–6× fructose-over-glucose HMF advantage. |
| `mp_holdout_glucose_alanine_130C_2h_pH50_Schibi` | HMF | **ANSWERED.** |
| `mp_holdout_glucose_alanine_130C_2h_pH80_Schibi` | HMF | **ANSWERED.** ⚠️ **The pH-5.0 and pH-8.0 arms will predict the SAME HMF**, because the module carries **no pH term on any HMF edge** — K5a declared gap G8: no paper in the cluster varies pH, so a pH term would have to be fitted across labs and matrices at once, which `k3` §B.2 forbids. The pH pair is therefore a **pre-registered structural miss**, and the size of the measured pH effect is the size of the gap. |
| both Schibilsky bundles | DMHF | **ANSWERED, and expected to UNDER-PREDICT BADLY — by one to three decades.** Reasoning, before scoring: the level is a declared transfer of a PENTOSE calibration onto a HEXOSE route, through a 1-deoxyglucosone pool that Martins' measured `k_odg_aa` drains ~10⁴× faster than the acetylformoin edge fills it. **If DMHF comes out within 3× of measurement, that is a coincidence and this wave will say so rather than claim it.** |

**DIRECTION, pre-registered for HMF: OVER-prediction.** The module has **no validated HMF
sink at cooking temperature** — Hamzalioglu's is clamped at 50 °C, self-degradation is a
single-temperature 0.9 %/7 d control carried with Ea = 0, and K5a G2 records that the
50–150 °C window is empty. A source-only node over-predicts.

**Aggregate commitment: at least 5 of the 7 rows become ANSWERED; the remaining refusals, if
any, must be SHARPER than the pre-B7 ones** (i.e. name a missing measurement rather than a
missing species).

---

## §4. WHAT WOULD FALSIFY THE BUILD

1. The fit RMS coming out **materially below 0.12 decades** — that would mean something was
   fitted beyond the one declared parameter (§2.2).
2. The two seeded starts disagreeing in `log10 k_dpo_af` by more than 1e-6 — the one free
   parameter is not identified.
3. `k_af_dmhf` moving the DMHF prediction by more than 1.2× over three decades — the
   "not rate-limiting" assumption is false and the constant is load-bearing (§2.3).
4. Any previously-answered exam row moving by more than 1.5× from the trunk perturbation
   (§1.3) — the furanic drains are larger than sized.
5. A grep finding any NaCl-arm literal, any coffee-arm literal, or any exam-side measured
   value inside `src/kinetic_core/*furanic*` or the frozen fit report.
6. `hmf_limb_shares` returning the same share for two systems with different sugar
   compositions — that would mean a branch fraction had been hard-coded after all.

---

## §5. WHAT THIS WAVE MAY NOT DO

* May not fit, tune or "estimate" **any activation energy for any furanone edge**
  (single-temperature literature, K5b §7.1).
* May not fit **anything** to Shu & Ho's 6.0 % GC area (K5b §8.6, by name).
* May not extrapolate **any Hamzalioglu constant above 50 °C** (K5a §7.3).
* May not ingest **any Gursul Aktag Table 2 activation energy** (0 of 43 reproduce; six are
  mathematically underivable) or **any Goncuoglu Tas activation energy** (0–1174 kJ/mol,
  six zeros, author-disclaimed twice, and the values are in a Table S4 that is not on disk).
* May not use **Lee et al. 2024's Arrhenius constants** — they are matrix→vapour transport
  and oven-wall deposition, not chemistry.
* May not hard-code **any branch fraction** between HMF's two limbs (K5a MUST-NOT #1).
* May not move any B1/B2/B3 fitted constant, any benchmark record, any hold-out role, or the
  declaration's existing text.
* May not run or trust a DFT barrier (owner policy, repository-wide).
