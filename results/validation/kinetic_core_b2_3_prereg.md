# Build Wave B2.3 — PRE-REGISTRATION

**Written and saved to disk BEFORE the refit ran, BEFORE the hold-out panel was
scored, and BEFORE the cutover exam was regenerated.** Every claim below is
scored in `kinetic_core_b2_3_fit_report.md` and
`kinetic_core_b2_3_holdout_report.md` whether it held or not. Claims that were
falsified are reported as falsified; the count of each is printed.

Mandate: `docs/reference/FIT_HOLDOUT_DECLARATION.md` **Amendment 9**, which
pre-declared this wave before it ran.

---

## 0. What this wave changed, in one paragraph

**No parameter was added, no rate constant, no fit row, no optimiser setting.**
The only structural change is a conservation fix. B2.2's post-freeze diagnosis
§3 found that sink steps delete a carboxylate- or amine-bearing solute without
saying where the titratable centre went, leaving the strong ion that titrated it
balancing nothing; Kang's 120 °C pot ran to a predicted pH 11.4 against a
measured 4.9. B2.3 carries the carboxyl into a new terminal pool (`CBX`),
declares the fate of the liberated amino nitrogen explicitly, and adds a fourth
**import-time** invariant — `ph_state.validate_charge_closure`, alongside
carbon, nitrogen and sulfur — that refuses to load the module unless every
step's titratable-centre movement is declared in the network's `CENTRE_LEDGER`
with its exact delta and a written basis.

### The audit's own finding, stated before the refit

Diagnosis §3 named **four** leaking sites (`r_cys_thermal`, `r_cys_h2s`,
`r_dpo_ddp`, `r_cys_actz`). Walking the whole network with the new audit finds
**fourteen**, of which the four largest by charge were *not* on that list:

| site | lane | centres deleted | named by §3? |
|---|---|---|---|
| `r_arp_dpo`, `r_arp_tdp`, `r_arp_dpo_thermal`, `r_arp_tdp_thermal` | sulfur | carboxyl + amine | **no** |
| `r_arp_decay` | sulfur | amine | **no** |
| `r_ttca_deg` | sulfur | amine (its carboxyl was in `ACID`, i.e. scaled by a fitted yield) | **no** |
| `r_cys_h2s`, `r_cys_thermal`, `r_dpo_ddp`, `r_cys_actz` | sulfur | carboxyl + amine | yes |
| `a_cys_sink`, `a_cys_glc` | acrylamide | carboxyl + amine | **no** |
| `r_fa_frag`, `r_aa_frag` | trunk | carboxyl (legitimate acid decomposition) | **no** |

The four Amadori enolisations are the dominant leak on Zhou's pot, whose entire
charge inventory is the 20 mmol/L of fed ARP they consume.

### The asymmetry, declared in advance

The carboxyl is **carried**; the liberated amino nitrogen is **declared
non-titratable at the point of release** (`ph_state.AMINE_FATE_BASIS`). This is
a chemistry claim, not a convenience, and it is the wave's strongest single
assumption. Its basis is Amendment 7's own declared FIT anchors: Zhou's pot
charges 40 mmol/L of amino nitrogen and finishes at a measured pH 3.22 / 3.42,
which 40 mmol/L of surviving ammonium (pK<sub>a</sub> 9.25) makes arithmetically
unreachable. It is also the textbook reason Maillard systems acidify at all.

**FULL DISCLOSURE OF THE ORDER OF OPERATIONS.** A pre-freeze probe
(`scratch/b23_encoding_probe.py`) was run at B2.2's *own frozen parameters*,
before this choice was made and before any refit, comparing the three admissible
encodings against Amendment 7's three FIT anchors and against Kang's 4.9 (an
**out-of-sample diagnostic**, declared in no column, already disclosed in B2.2's
diagnosis §3). It returned:

| encoding | Zhou 6/7/8 cooled pH (anchors 3.22 / 3.42 / 5.07) | Kang 120 °C (diagnostic 4.9) |
|---|---|---|
| B2.2 as shipped | 3.29 / 3.49 / 5.07 | **11.57** |
| carry carboxyl, declare amine destroyed — **SHIPPED** | 2.94 / 3.02 / 3.47 | **4.85** |
| carry both centres | 5.04 / 5.11 / 5.89 | 9.11 |

The choice was made on the arithmetic of the FIT anchors and on the acidification
argument. That it also lands Kang's diagnostic at 4.85 against 4.9 arrives as
corroboration, and this paragraph exists so that it cannot later be presented as
the reason.

---

## 1. THE FIT (F-block)

Scored against `kinetic_core_b2_3_fit_report.json`. B2.2's frozen values are
χ²ᵣ **1.80** on 58 rows / 48 parameters, `acid_yield` **0.968**,
`arp_secondary_ammonium_pKa` **10.77**, `Ea_decay_thiol_sink` **248.0** (its
ceiling is 250), `Ea_decay_carbonyl_sink` **150.2**.

| # | claim |
|---|---|
| **F-1** | χ²ᵣ lands between **1.5 and 4.0**, and is **worse** (larger) than B2.2's 1.80. Reason: the fix removes the optimiser's cheapest escape from the leak without giving it anything back, and the fit row count and parameter count are both unchanged. |
| **F-2** | The three Zhou drift endpoints are still reproduced within **0.5 pH units** each. This is the fix's own calibration target and it should survive; if it does not, the amine declaration is wrong and the wave says so. |
| **F-3** | `arp_secondary_ammonium_pKa` falls **below 9.5**, i.e. away from B2.2's 10.77 and toward the 7–8 an N-glycosylated secondary amine should give. Reason: B2.2's diagnosis §5 contradiction 2 argued that 10.77 is a nuisance value the fit chose in order to make the Amadori electrically invisible — which is precisely how an optimiser stops a leak it cannot repair. With the leak repaired, the incentive is gone. **This is the wave's sharpest falsifiable prediction: if the pKa stays above 10, the contradiction-2 diagnosis was wrong.** |
| **F-4** | `acid_yield` falls **below 0.90**. Reason: `CBX` now supplies full-strength acid that `ACID` used to have to supply through a fitted fraction. |
| **F-5** | `Ea_decay_thiol_sink` stays **above 200 kJ/mol**, i.e. still bound-limited. The mechanism B2.2's diagnosis §3 identified (a sink whose *absolute* rate must be tiny across 100–145 °C) is untouched by a charge fix. |
| **F-6** | The number of individually identified parameters stays in the range **2–8 of 48**. The fix adds no information to the objective. |
| **F-7** | The Kang 100/120 °C pots no longer run alkaline: both **cooled endpoints land below pH 7.0**, against B2.2's 11.42. (Kang's 4.9 is a diagnostic, not a scored row; the claim is about the *direction*, and it is a claim the conservation law makes on its own.) |
| **F-8** | The Hofmann 0.5 M phosphate pot still drifts — its van Slyke capacity at pH 5 is ≈15.7 mmol L⁻¹ per pH unit against a modelled acid load in the tens — so **`hofmann_pentose_pH5` finishes below pH 4.5 cooled**. Buffers in this corpus do not hold, and the fix does not change that arithmetic. |

---

## 2. BLIND RE-SIT (a) — the full B2.1/B2.2 hold-out panel

Trajectory to date: **B2.1 15/33 → B2.2 12/33**. Amendment 9 demotes
`kang_140C_MFT` to seen-diagnostic **permanently**, so B2.3's own denominator is
32; both totals are reported and the like-for-like is printed first.

| # | claim |
|---|---|
| **H-1** | Gating lands between **10 and 16 of 33** on B2.2's denominator. Wide, and deliberately so: a charge fix moves the operating pH of every unbuffered pot, and every row scored against a nominal pH label is scored at a different operating point. |
| **H-2** | Between **3 and 12 rows change status** in either direction. B2.2 moved 6 FIXED / 9 REGRESSED on a smaller structural change. |
| **H-3** | The **Zhou pH-8 block improves**, in aggregate fold error, relative to B2.2. Reason: the pH-8 pot has the largest SID and therefore the largest bookkeeping leak, so it is the row block the fix should touch most. |
| **H-4** | The **`kang_switch_on_*` rows still fail**. The switch-on is a barrier-ratio phenomenon; nothing in a charge fix creates one. |
| **H-5** | `cerny2007_branch_responsiveness` still fails. B2.1 proved its ceiling is arithmetic. |
| **H-6** | `zhou_dimer_share_pH_invariance` and the three `zhou_120C_dimer_share_pH*` rows are **not** collectively resolved — at most 2 of the 4 pass. The Kumazawa-vs-Zhou tension B2.1 and B2.2 both recorded as unadjudicated is not a charge-conservation problem. |
| **H-7** | `vanseeventer_50C_MFT_per_day` still passes and `hofmann2002_brew_80C_FFT` still passes. Both are low-temperature rows in systems with little acid production, so the fix should barely reach them. |

---

## 3. BLIND RE-SIT (b) — the cutover exam, BOTH WAYS

Amendment 9 clause 2 requires the exam to be reported **both** with the completed
buffer fields and **as-was**, permanently, in one artifact. B2.2's exam numbers:
core paired median **10.65×**, old lane **12.65×**, Hofmann 145 °C family median
**4.16×** (4/10), Yiltirak **0/8 median 262.3×, worst 1475×**, 17 declensions.

| # | claim |
|---|---|
| **X-1** | The **buffer-completed** Yiltirak family median is **at least 3× better** than the as-was median. B2.2's §4 sized the buffer field alone at an 11× swing in predicted FFT on exactly this system shape. |
| **X-2** | Yiltirak nonetheless stays **0/8** buffer-completed. An 11× correction against a 262× median does not close it, and the residual is a chemistry failure that the schema gap was masking. |
| **X-3** | The **as-was** column is *not* uniformly worse. Some rows will be closer as-was, by cancellation of two errors, and reporting only the completed column would hide that. |
| **X-4** | The core's paired median stays in the band **5×–25×** in both columns. |
| **X-5** | The core **still beats the old lane** on the paired median in the buffer-completed column. |
| **X-6** | The 17 declensions stay 17. Nothing in this wave adds or removes a declension. |
| **X-7** | The Chang two-arm defect is **not** fixed: both arms still predict the same value. The acrylamide lane has no pH term and B2.3 does not give it one. |

---

## 4. BLIND RE-SIT (c) — Kang 140 °C, scored but NOT COUNTED

Per Amendment 9's trigger disclosure, `kang_140C_MFT` is **permanently**
seen-diagnostic and its B2.2 nominal pass is ruled not-counted. It is still
scored and printed.

| # | claim |
|---|---|
| **K-1** | The Kang 140 °C **pH trajectory is no longer six units wrong**: the cooled endpoint lands below 7.0. |
| **K-2** | `kang_140C_MFT`'s fold error **degrades** relative to B2.2's 2.15×, to somewhere between 2× and 30×. B2.2's diagnosis §3 argued that its pass was a coincidence of a broken trajectory; repairing the trajectory should therefore cost it. **This is the second sharpest prediction in this file, and it predicts the wave's own number getting worse.** |
| **K-3** | `kang_140C_FFT` stays failing. |
| **K-4** | `kang_140C_cys_conversion` stays passing — it is governed by Kang's own immovable 55.1 kJ/mol. |

---

## 5. BUFFER-FIELD COMPLETION (Task 2) — claims about the DATA, not the model

| # | claim |
|---|---|
| **B-1** | At least **6 of the 21** bundles resolve to `buffer_unknown`, because their source paper is not on disk. |
| **B-2** | At least **one** bundle's paper affirmatively states **no buffer** (plain water), which is a positive finding and must be recorded as `buffer: none`, not as unknown. |
| **B-3** | The `measured_volatiles` block of **every one of the 21 bundles is byte-identical** before and after the edit, proven by a SHA-256 hash test committed alongside. Compound lists, `evidence_class` and roles likewise. |
| **B-4** | At least one bundle's recorded `conditions.ph` **disagrees with its source paper**, and the disagreement is reported rather than silently corrected — editing an executable condition would change a prediction. |

---

## 6. What would make this wave a failure

Stated in advance so that it cannot be renegotiated afterwards:

1. **F-2 falsified** — the drift anchors no longer reproduce. That would mean the
   amine declaration is wrong, and the wave would have to report that its central
   claim does not survive its own calibration target.
2. **H-1 below 10/33.** A conservation fix that costs three or more net gating
   rows is a fix whose accompanying assumption is doing damage, and the honest
   response is to report it, not to tune.
3. **K-1 falsified** — the Kang pot still runs alkaline. That would mean the leak
   was not the cause, and B2.2's diagnosis §3 was wrong about the mechanism.

None of these would be repaired by changing a parameter after seeing the score.
If any occurs it is written up in `kinetic_core_b2_3_diagnosis.md` as a post-hoc
explanation, filed as a post-hoc explanation.
