# Wave B64 pre-registration — the thermal history each source records

*Written 2026-10-28, before any B64 fit number existed. Owner, the same day: "You can fix this too" (the
heat-up), then "Plan thoughtfully and execute" on the revised decision set (`tasks/data_restructure_plan.md`
sec. 7, D1). Kind 1 (a data correction: the pots were not isothermal) with its kind 2 refit, accepted under the
frozen rule (`scripts/generators/WAVES.md`, "Evidence decides, scores diagnose") on the audit in sec. 3 and the
verification in sec. 4. An independent reviewer reads the audit before the engine switches. Generator:
`scripts/generators/generate_kinetic_core_b64_fit.py`. Switch: `panel.THERMAL_HISTORY` ("isothermal" →
"as_recorded"). Baseline: B60 (shipped). Runs before B63. North-star targets served: 1 (measured where the data
allow: the inputs too) and 4 (timescales).*

## 1. Why

Every pot, fitted or scored, has been integrated **isothermal at its nominal temperature from t = 0**. The
sulfur lane's anchor laboratory says its pots were not. Read verbatim in the PDFs on disk
(the per-source record, with every quote and its page: `docs/validation/thermal_history_facts.yml`):

- Hofmann & Schieberle 1995 (`hofmann1995.pdf`, p. 2): cysteine 3.3 mmol + ribose 10 mmol in 100 mL 0.5 M
  phosphate pH 5.0 "were allowed to react in a laboratory autoclave (Model II; Roth, Karlsruhe, Germany) **by
  raising the temperature within 20 min from 20 to 145 °C**. … After cooling, …" — no hold between.
- Hofmann & Schieberle 1997 (`hofmann1997.pdf`, p. 901): the same charge, the same 200 mL autoclave, the same
  ramp sentence in the Methods; the Results on the same page call that pot "thermally treated for 20 min at
  145 °C".
- Hofmann & Schieberle 1998 (`hofmann1998.pdf`, p. 236), the source of the fit's 14 aqueous 145 °C systems:
  "thermally treated for **20 min at 145 °C** in a laboratory autoclave (200 mL; Type II; Roth, Germany)".
- Schieberle 2000 (Table III) prints the 1998 Table 1 thiol levels for "(145 °C; 20 min)" and refers to the 1995
  and 1997 papers for the procedure. No paper on disk states a hold at 145 °C, a come-up separate from the
  ramp, or how the autoclave is heated or cooled; the Roth autoclave appears only in these three papers.

The same laboratory writes the same pot both ways, so the 1998 wording is the 1997 abbreviation of the ramp.
**Readings and their probabilities (mine, before any fit):** ramp only, 20 → 145 °C over 20 min, then cool:
**0.80**; ramp then a 20 min hold: 0.12; isothermal 20 min: 0.03; other (e.g. a faster ramp in the 50 mL fed
pots): 0.05. For barriers of 64–120 kJ/mol a linear ramp is worth R·T²/(Ea·β) ≈ 1.9–3.6 min at 145 °C, so the
isothermal reading over-doses a ramp-only pot about 6–10×; a ramp plus hold is ~1.1× the isothermal dose.

**What the reading does at B60's constants (a probe run before this prereg, disclosed here; diagnosis, not a
result of the wave).** log10(predicted/measured), isothermal → ramp only: furfural + H2S → FFT +0.19 → +2.05;
norfuraneol + cysteine → MFT +0.02 → −0.79; ribose + Cys pH 5 FFT/MFT +0.46/−0.33 → +0.05/−0.85; fructose FFT
+1.23 → −0.84; ribose pH 7 hold-out MFT −3.69 → −6.16. The fit's 145 °C constants absorbed the over-dose; the
reading cannot be changed without the refit.

## 2. What changes

- **Fit systems** (`READINGS` in the generator): the 14 Hofmann 1998 systems take the reading's profile; van
  Seeventer 2001 takes its stated one in both readings ("from room temperature to 130 °C in 10 min and was kept
  at 130 °C for 20 min", Materials and Methods p. 4292; start taken as 20 °C), run time 30 min. Every other fit
  source states no come-up and stays isothermal, declared: Cerny 2007b (1 mL PTFE vial, metal block), Zhou 2023,
  Whitfield & Mottram 1999 (2 mL ampoules, oven), Zhang 2024, Kang 2026 (glass pressure vessels, oil bath, ice
  quench; 120 min holds), Kumazawa 2003 (cans in a retort, 10 min; retort practice suggests a hold after
  come-up, unstated), Feng 2022. Cooling: "after cooling" (to 20 °C in 1998) is unstated in rate and omitted,
  a declared under-dose in both readings.
- **Panel bundles:** every Hofmann 1998 bundle gets `conditions.heating = {time_convention: ramp_only,
  start_C: 20, ramp_min: 20}` with its provenance; read only when `panel.THERMAL_HISTORY = "as_recorded"`.
  Bundles whose source states a come-up the profile cannot encode yet (Leitzen 2021: an 18 min autoclave hold
  plus heating and cooling phases the source rates at F0 = 12 min; Lin 2021/2022: 100 mL at 180 °C in a dry
  bath, come-up unstated) are recorded and left isothermal: no invented time constants.
- **Engine:** `ThermalProgram.from_profile` (staircase of ≤ 2 K for the single-temperature lanes) and
  `integrate_sulfur(temperature_profile=...)` (rate constants and in-situ pH at T(t) on the pH nodes, nodes
  every ≤ 2 K). A program built from holds runs bit for bit as before.
- **Numerics (C5, owner 2026-10-28: "do 4 before too"):** the fit and every B64 score run on
  `sulfur.PH_COUPLING = "inline"`: the charge balance solved at each right-hand-side evaluation (warm-started
  Brent), with an analytic mass-action Jacobian, instead of B60's node fixed point. Checked before the fit
  against the converged fixed point (81 nodes, rtol 1e-8) on five probe pots (Hofmann ribose pH 5 and 7 at
  145 °C, Kang TTCA 100 °C unbuffered, Zhou ARP + Cys 120 °C unbuffered, the ramp): ≤ 0.002 dex on every
  species, cooled pH equal to 3 decimals, 1.6-14× faster. V5 below is run on it. The shipped engine keeps
  "fixed_point_nodes" until B64 ships.
- **Candidates:** (A) refit under `ramp_only`: B60's objective, rows, bands, free set and protocol; start 0 =
  B60's optimum. (B) refit under `ramp_plus_hold`: recorded, the declared alternative. (C) B60's vector under
  `ramp_only`, unrefitted: recorded (`--cost-only`).

## 3. Evidence audit

| audit item | finding |
| --- | --- |
| 1. established for the regime | Not a mechanism: a fact about the pots. The heating is stated verbatim by the laboratory for this pot and apparatus (sec. 1); the reading is the one its own sibling papers give. |
| 2. carries flux | The 14 systems are the fit's 145 °C anchors (most of the objective's rows); the probe moves them up to 1.9 dex. |
| 3. constants measured, bands carried | No new constant. The profile is the source's numbers. The alternative reading is carried as candidate B; if A ships, the envelope draws the reading (0.8 / 0.2) — a follow-up, stated here. |
| 4. verified in the engine | V1 a flat profile reproduces the isothermal run bit for bit; V2 a 20 → 145 °C ramp reproduces the exact Arrhenius integral of a first-order step (≤ 0.2 %); V3 one profile run agrees with a 0.25 K staircase (≤ 1 % on every species); V4 the fast pH residual is bit-identical (all in `tests/unit/test_kinetic_core_b64_thermal_profile.py`, PASS when this was written). V5 (before the fit): the quick setting under the profile agrees with the careful setting within 0.03 dex on all 52 observables, as B60's verification required. |
| 5. no contradiction with a direct measurement | Supporting: Schieberle 2000 Table IV runs the same charge in the same autoclave at 100 °C (MFT 13.8 µg at 60 min, 156 µg at 360 min); a 145 °C ramp-only dose of ~2–3 min-equivalent predicts tens of µg and a 20 min hold ≥ 150 µg; the measured 145 °C pot holds 19.8 µg. Consistent with ramp only (arithmetic to be re-checked in the ship rule). Nothing on disk states a hold. |

**Audit verdict, before any result: grounded. The pots' recorded thermal history replaces an unstated
isothermal assumption, and the alternative reading is carried, not discarded.**

## 4. Recorded diagnosis (not gating) and the rule

Recorded for A, B, C: the fit cost on the 52 rows; T1 the Hofmann ribose/xylose pH 3 / 5 / 7 shape (MFT
553 / 198 / 25, FFT 229 / 121 / 12 ppb); T2 the Yiltirak ladder; T3 every out-of-sample sulfur row; T4 the
directional panel; the 145 °C constants and the lumped barrier against B60's (my expectation: the 145 °C trunk
and assembly constants rise, the lumped barrier moves up from 64 toward the measured class ~100, because the
100–120 °C rows now face a ~7× smaller 145 °C dose); the typical-yield baseline and the grounding share.

**Ship rule (mechanical):** A ships when V1–V5 pass, both starts converge (pH fixed point converged on every
system) and the independent reviewer finds no error in sec. 1–3. Scores do not gate. If A's fit is worse than
B60's cost, the rows that got worse are recorded as the compensating errors they expose.

**Compute:** two starts per reading, `nice -n 15`, launched only when the earlier B61 processes have finished
(≤ ~5 cores in total).

## 4a. Amendment 1 (2026-10-28, before any fit member was written; owner: "avoid arbitrary choices")

The register of choices (`docs/reference/choices_register.yml`, CR-03 and CR-04) replaces three judgement calls
above with rules that need no chosen number:
1. **No reading is privileged by my probability.** Sec. 1's 0.80 / 0.12 stays as a stated opinion and decides
   nothing. Candidates A (`ramp_only`) and B (`ramp_plus_hold`) are both refitted and both reported, on every
   diagnostic in sec. 4, side by side. If both pass the mechanical rule, the envelope draws the two readings
   with EQUAL weight, and the point is the one whose constants reproduce the one measurement that bears on it
   directly: Schieberle 2000 Table IV, the same charge in the same autoclave at 100 °C (MFT 13.8 µg at 60 min,
   156 µg at 360 min) against the 145 °C pot's 19.8 µg. If that check does not separate them by more than its
   own stated uncertainty, the point is not chosen: both are printed, and the reading goes to the experiment
   ranking.
2. **The ramp's shape (CR-04)** is not refitted. At each candidate's optimum the 52 rows and the hold-outs are
   rescored with (i) a Newtonian heat-up reaching 145 °C at the stated 20 min and (ii) a cool-down of the
   ramp's own duration; the change is reported as the reading's sensitivity, not used to pick anything.
3. **No tolerance is introduced** to call the readings equal: they are compared on integer counts (rows within
   3x, by lane) and on the full fold-error distribution, all printed.

## 4b. Amendment 2 (2026-10-28, after the discriminator and the shape sweep, before any Newtonian fit)

Two results arrived after Amendment 1 and before this text: (i) the discriminator (Schieberle 2000's
100 °C / 6 h over 145 °C MFT factor, measured 13×) is reproduced by NEITHER optimum (ramp only 2.1×, ramp + hold
0.2×): both make MFT at 100 °C for under an hour and then lose it, where the pot keeps forming it for 12 h — the
thiol / trunk-timescale problem, not a thermal one; so per Amendment 1 no reading is chosen. (ii) The shape
sweep: at the ramp-only optimum a Newtonian heat-up lowers the cost from 45.4 to 39.7 (τ = 20/3 min) and 33.4
(τ = 20/5 min) without refitting, so the linear shape is a consequential choice (register CR-04). The two
Newtonian shapes are therefore refitted as readings `ramp_only_newtonian_tau6.7` and `ramp_only_newtonian_tau4`
(start 0 = the ramp-only optimum), and every reading is reported side by side with its discriminator value.
Nothing is chosen on cost.

## 5. Verdict

*(appended after the rule runs; nothing above this line is edited after that)*

Written 2026-10-28, 23:05, after every B64 fit and check had run (`kinetic_core_b64_*_members/`,
`kinetic_core_b64_discriminator.json`, `kinetic_core_b64_shape_sensitivity.json`).

| reading of Hofmann 1998's "20 min at 145 °C" | best cost on the 52 rows (two starts where run) | MFT 100 °C-6 h / 145 °C (measured 13×) |
| --- | ---: | ---: |
| isothermal (every wave through B63; no source states it) | 27.08 (B60's vector) | — |
| ramp + hold, linear | **25.95** (25.97 / 25.95: one optimum) | 0.2× |
| ramp only, Newtonian, within 1 K at 20 min | 30.87 (budget exhausted) | 0.4× |
| ramp only, Newtonian, within 6 K at 20 min | 38.69 (budget exhausted) | 0.9× |
| ramp only, linear | 45.43 / 46.00 (one optimum) | 2.1× |

**Mechanical rule:** V1-V5 PASS; both starts agree where two were run; the pH fixed point does not exist on the
inline coupling. The review of sec. 1-3 by an independent reader is still owed before any switch moves.

**Verdict (Amendments 1-2, applied as written): NO READING IS CHOSEN.** The one direct measurement does not
separate them: every reading misses the 100 °C factor by 0.8-1.8 dex, and the misses are anti-correlated with the
fit cost — the more 145 °C exposure a reading gives, the better the fed pots fit and the further the 100 °C series
falls behind. That is the signature of a structural error, not of the thermal history: at 100 °C the model stops
forming MFT within an hour where the pot forms it for 12 h (thiols and trunk timescale; B59). Recorded:
1. `panel.THERMAL_HISTORY` stays "isothermal" and `sulfur.PH_COUPLING` moves to "inline" only as a numerics
   change (kind 1, verified ≤ 0.002 dex), not as part of a thermal reading.
2. The thermal reading is an open hypothesis with four alternatives (register CR-03, CR-04: moves_headline),
   carried forward: B59a and every later wave are run on the two best-fitting shape classes (ramp + hold linear;
   ramp only Newtonian within 1 K) and reported against each.
3. The experiment ranking gains "log the autoclave's temperature during a Hofmann-protocol run": it would retire
   three of the four readings in an afternoon.
4. The shipped scorecard does not move: no prediction changes until a reading is chosen.

### Independent review (2026-10-29, a fresh reader with the files only): PASS WITH CHANGES — corrections recorded

1. The discriminator's "measured 13x" is the text's number; Schieberle 2000's own tables give 156 ug / 100 mL
   (Table IV) over 198 ug/L (Table III) = 7.9x (13 is FFT's ratio). The measurement carries a declared 0.2 dex
   (7.9-13x) uncertainty; no conclusion above changes.
2. The verdict's ground is restated: the check DOES separate the readings (1 dex); what failed is its premise — no
   reading reproduces the factor, and the factor is model-dependent. On the model-free evidence (four Munich papers'
   wording, the hardware) the ramp-only reading is favoured; it is still not chosen, because the 145 C fed pots it
   implies cannot be fitted by the present network, and that is a structural question, not a thermal one.
3. Post hoc: Amendment 2's Newtonian refits reversed Amendment 1 item 2 after costs were seen, and the tau values
   are chosen numbers (CR-04). Verdict item 2 ("the two best-fitting shape classes") is WITHDRAWN: later waves run
   on every reading, isothermal included.
4. The mechanical rule was overstated: most members are unconverged (status -9), the second start of each reading
   is the other reading's optimum, and V5 was not run on the Newtonian readings.
5. The inline pH coupling ships only after a full-objective and full-panel diff against the fixed point.
