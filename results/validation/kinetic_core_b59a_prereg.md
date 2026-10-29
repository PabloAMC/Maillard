# Wave B59a pre-registration — the Strecker sulfur released as H2S

*Written 2026-10-28, before any number with the switch on existed. Approved by Pablo as the first half of B59
("Go ahead with B59 after B61 and B62"; reordered to follow B64 on 2026-10-28, decision D2: on the literature
alone). Accepted under the frozen rule (`scripts/generators/WAVES.md`, "Evidence decides, scores diagnose") on
the audit in sec. 3. Switch: `parameters_sulfur.STRECKER_SULFUR`. Reaction: `r_dpo_ddp_h2s`. Generators to
come: `generate_kinetic_core_b59a_fit.py`, `..._ship_rule.py`. North-star targets: 1 (a route the literature
draws, at a measured share) and 4 (it changes what the pots can do over time, not only branch ratios).*

## 1. Why

`r_dpo_ddp` (1-deoxypentosone + cysteine, Strecker) sends cysteine's sulfur, as mercaptoacetaldehyde, to the
inert `FRAG_S`. On B60 it consumes 98-100 % of the cysteine in the Hofmann ribose pots in seconds
(`SULFUR_FLUX_MAP.md` sec. 2), so at pH 7 the pot makes 0.0002 mM H2S and MFT is ~4 900× low while 13.5 mM
norfuraneol waits for sulfide. The textbook route — mercaptoacetaldehyde → H2S + acetaldehyde — is absent.

## 2. What changes

- A second channel with the same reactants, `r_dpo_ddp_h2s`: {DPO, Cys} → {DDP, CBX, 2 FRAG_C, FRAG_N, H2S}.
  Its constant is `h2s_share × k_dpo_ddp`, and `r_dpo_ddp` keeps `(1 − h2s_share) × k_dpo_ddp`, so the Strecker
  step's total rate, its barrier and its pH law are unchanged; only where the sulfur goes moves. Off: the
  pre-B59a network exactly (the new rate is 0).
- **The share is not fitted and not chosen.** Kobayasi & Fujimaki 1965 measure H2S as 4.2 % and 8.0 % of the
  cysteine CHARGED (100 °C, 1 h, 1 or 2 equivalents of glucosamine as the reductant). A yield per cysteine
  charged is a LOWER BOUND on the share of Strecker-processed sulfur released. The wave runs at both measured
  values (alternatives `h2s_share` = 0.042 and 0.080; register CR-11) and reports both. A higher share is not
  excluded by anything on disk; that is stated, not hidden.
- **Baseline:** B64 does not choose a thermal reading (its Amendment 1). B59a is therefore run on the two B64
  readings with the best fit per shape class — `ramp_plus_hold` (linear) and the best Newtonian ramp-only — so
  the result is reported against each, not against a reading chosen for it.
- **Candidates per baseline:** (A) refit the shipped free set with the switch on (kind 2); (B) the baseline's
  vector with the switch on (recorded).

## 3. Evidence audit

| item | finding |
| --- | --- |
| 1. established for the regime | Cysteine + carbonyl (ninhydrin, glucosamine) boiled in water at 100 °C gives mercaptoacetaldehyde, H2S and acetaldehyde (Kobayasi 1965, p. 698); the H2S is released by the enaminol the Strecker step itself makes. Du 2023's labelling (glucose + cysteine + aldehydes) attributes the H2S and NH3 of its thiophenes and pyridines to cysteine's own Strecker degradation. Same catalysis (amine-carbonyl Strecker), aqueous, 100-145 °C. |
| 2. carries flux | `r_dpo_ddp` is the main cysteine sink in every reference pot on B60 (98-100 % at Hofmann pH 5 and 7, Zhang 2025b). |
| 3. constants measured, bands carried | No new fitted constant. The share is the source's two measured yields, carried as alternatives; the Strecker constant keeps its existing fit. Declared transfer: glucosamine (an enaminol-forming amino sugar) and ninhydrin stand in for the 1-deoxypentosone; risk direction unknown. |
| 4. verified in the engine | V1: switch off is bit-identical on the full objective. V2: with the switch on, the two channels' constants sum to the old one at every T and pH, and the sulfur balance closes. V3: the engine reproduces Kobayasi's H2S share in its own description of a cysteine + DPO pot only in the sense the source licenses — as a lower bound (H2S / cysteine charged ≤ the share). |
| 5. no contradiction with a direct measurement | Mizutani 1960 (H2S and NH3 from cysteine thiazolidines at 120 °C) and Edwards & Wedzicha 1992 (cysteine + 3-deoxyhexosulose: fast thiazolidine, slow release) bear on cysteine's other fates, not on this share; neither contradicts it. Zheng & Ho's cysteine thermolysis (already in the model) stays as the second H2S source. |

**Audit verdict, before any result: grounded as a route with a measured lower bound on its share.**

## 4. Recorded diagnosis (not gating)

The Hofmann ribose pH 3 / 5 / 7 shape (MFT 553 / 198 / 25, FFT 229 / 121 / 12 ppb; the pH 3 and 7 arms are
hold-outs); H2S made at pH 7; the Schieberle 2000 100 °C series and its 100 °C-6 h / 145 °C MFT factor (13×);
every out-of-sample sulfur row; the contrast diagnostic; the typical-yield baseline; the share of reference-pot
flux on measured constants. Expectation (mine): pH 7 MFT rises by orders of magnitude, pH 3 changes little,
FFT rises with the extra H2S, and the 100 °C factor moves toward 13 only if the cysteine lasts longer, which is
B59b's job, not this one's.

## 5. Verdict

*(appended after the rule runs; nothing above this line is edited after that)*

Written 2026-10-28, 23:50, after the four refits (2 baselines × 2 shares; 1 200 evaluations on ramp + hold, 600
on the Newtonian baseline) and the diagnosis (`kinetic_core_b59a_diagnosis_{ramp_plus_hold,ramp_only_newtonian_tau4}.json`).

| baseline (B64 reading) | arm | cost (52 rows) | MFT log10 pred/meas pH 3 / 5 / 7 | FFT pH 3 / 5 / 7 | H2S made, pH 7 |
| --- | --- | ---: | --- | --- | ---: |
| ramp + hold | switch off | 25.95 | −2.62 / +0.05 / **−6.85** | +0.53 / +0.97 / −5.31 | 2e-8 mM |
| ramp + hold | share 0.042 | 32.36 | −2.94 / −0.06 / **+0.32** | +1.37 / +1.73 / +1.03 | 0.46 mM |
| ramp + hold | share 0.080 | 34.78 | −3.13 / −0.03 / **+0.51** | +1.83 / +1.90 / +1.07 | 0.89 mM |
| Newtonian ramp only (τ = 4 min) | switch off | 30.87 | −2.81 / −0.24 / **−6.02** | +0.09 / +0.52 / −4.30 | 1e-7 mM |
| Newtonian ramp only | share 0.042 | 41.59 | −3.16 / +0.10 / **+0.51** | +0.36 / +1.49 / +1.67 | 0.82 mM |
| Newtonian ramp only | share 0.080 | 43.95 | −3.30 / +0.20 / **+0.68** | +0.60 / +1.68 / +1.72 | 1.58 mM |

**Audit (sec. 3): stands.** V1-V2 PASS (switch off: the new channel is exactly zero; on: the two channels sum to
the old constant, C/N/S close to 1e-9). V3 holds as a bound. No direct measurement contradicts the route.

**Diagnosis.** On both baselines and both measured shares the never-fitted pH 7 arm moves from 6-7 decades low to
within 2-5× (MFT), and pH 5 stays within 2×: the pH 7 failure was H2S starvation, as the flux map said. The fit cost
rises because the rows now expose the errors the starvation hid: MFT at pH 3 ~3 decades low (the enolisation pH
law; B59b), FFT 2-50× high (the furfural branch; B62), and the 100 °C / 145 °C factor (0.09-0.10 against 13; the
cysteine lifetime; B59b). Those are recorded as the compensating errors B59a exposes, not as reasons against it.

**Verdict: SHIP, pending the independent review**, with the share unresolved: both measured lower bounds are
carried (register CR-11), the point at the lower bound (0.042, the one equivalent of reductant Kobayasi measured),
the envelope across both. The baseline stays unresolved per B64: B59a ships on top of whichever thermal reading
B64 ships with, and both baselines' constants are recorded here.

### Independent review (2026-10-29, a fresh reader with the files only): FAIL — the verdict above is WITHDRAWN

Accepted in full: (1) Kobayasi's H2S per cysteine charged is not a lower bound on a Strecker share (no blank,
glucosamine is a reductant not a deoxyosone, cysteine itself may be cleaved, and H2S scales ~1.9x with the
reductant charge): a fixed share is the wrong functional form; the release is reductant-limited. (2) REDOX ERROR in
`r_dpo_ddp_h2s`: MAA + 2H -> acetaldehyde + H2S takes the two electrons from the enaminol, so the C5 fragment must
leave at the deoxyosone oxidation level, not as the 2e-reduced DDP the channel emitted; C, N, S balanced, H and
redox did not, and DDP feeds the direct MFT route, so the pH 7 "rescue" above is not evidence until this is fixed.
(3) Du 2023 traces carbon, not sulfur; Wedzicha 1991 (the released thiol persists 25 h at 55 C) was omitted. (4) V3
was vacuous; no mechanical ship rule was pre-registered; fits unconverged; no isothermal arm (its refits were
stopped on this review). The switch stays OFF. Redesign before any rerun: a reductant-limited, redox-balanced step,
its constant treated as an order-of-magnitude transfer with unknown direction of risk, a real engine check, a
pre-registered ship rule, converged fits on every B64 reading including isothermal.

## B59a-r — revised pre-registration (2026-10-29, after the review and before any number with the revised channel)

**Change.** `r_dpo_ddp_h2s` now returns the 1-deoxyosone: {DPO, Cys} → {DPO, CBX, 2 FRAG_C, FRAG_N, H2S}. Chemistry:
the Strecker step gives mercaptoacetaldehyde + CO2 + the aminoreductone (the dicarbonyl having taken the amine and
2 H); the aminoreductone reduces the mercaptoacetaldehyde to H2S + acetaldehyde and, re-oxidised to the imine,
hydrolyses to the same deoxyosone + NH3. Net Cys + H2O → H2S + CH3CHO + CO2 + NH3; H, O and electrons balance. The
release is reductant-limited by construction (≤ one H2S per Strecker event).

**The share.** The fraction of aminoreductones that reduce mercaptoacetaldehyde rather than going on to DDP is not
measured. Kobayasi's 4.2 / 8.0 % H2S per cysteine charged (100 °C, 1 h, 1 / 2 equivalents of glucosamine as the
reductant) is carried as an ORDER-OF-MAGNITUDE TRANSFER with unknown direction of risk — not a bound (the review:
no blank; glucosamine is not a deoxyosone; cysteine itself may be cleaved). Both values are run (CR-11); a value is
chosen only if a measurement in the fit rows identifies it (profile), otherwise both are reported.

**Evidence audit, revised.** 1: Kobayasi 1965 (cysteine + carbonyls at 100 °C give mercaptoacetaldehyde, H2S and
acetaldehyde; reductive mechanism); Du 2023 is withdrawn as evidence for the sulfur (it traces carbon). 2: unchanged
(the Strecker step is the main cysteine sink). 3: no new fitted constant; the share is a declared transfer.
4 (V3, now a reproduction): in the engine, DPO 10 mM + Cys 10 mM, unbuffered pH 6 start, 100 °C, 1 h — the H2S made
per cysteine charged is recorded beside Kobayasi's 4.2-8.0 % (same order of magnitude required at either share; a
check of the transfer, not a fit). 5: Wedzicha & Edwards 1991 (cysteine + 3-deoxyosone at 55 °C: the released thiol,
proposed as mercaptoacetaldehyde, persists ~25 h) is a direct measurement against fast H2S release at 55 °C; the
step's barrier (the Strecker step's) makes it ~10^2-10^3 slower at 55 than at 100-145 °C, so no contradiction at
cooking temperature — declared, not resolved.

**Mechanical ship rule (written now).** Ship when ALL hold: (a) V1-V3 pass; (b) an independent reviewer finds no error
in this section; (c) every fit is run on EVERY B64 reading (isothermal, ramp + hold, ramp only linear, Newtonian
τ = 4 and 6.7 min) at both shares, each to the optimiser's own convergence or to 1 200 evaluations, whichever first,
and reported as converged or not. Scores do not gate; the shipped engine takes the isothermal-reading constants at the
share the profile identifies, or both shares as alternatives if it identifies neither.

**V3 (run 2026-10-29, at B60's constants, before any B59a-r fit):** DPO 10 mM + Cys 10 mM, unbuffered from pH 6,
100 °C, 1 h → H2S 0.69 % (share 0.042) / 1.23 % (share 0.080) of the cysteine charged, against Kobayasi's 4.2-8.0 %:
the same order of magnitude, low by 3.4-6×. PASS as written (order of magnitude), recorded as low.

## B59a-r — Amendment 1 (2026-10-29, before any corrected run)

The eight B59a-r refits and the two ramp-only ones (costs 34.5-56.8 against 26-27 for the baselines; 9 of 10
budget-exhausted) freed all 23 coordinates. The frozen rule refits only the coordinates the data identify
(`WAVES.md`). At the isothermal optimum, share 0.042, the thiol-sink barrier (identified, but loosely: σ 33 kJ/mol) fell from 101 to
59 kJ/mol with `k_thiol_decay` +0.85 dex; the Schieberle 100 °C factor fell to 0.01 (measured 7.9-13) and the
Nestlé shape inverted (MFT at 30 min 5.6× its 4 h value). These members are withdrawn to `withdrawn_all_free_set/`
as a recorded diagnosis, not as the wave's result. B59a-r is re-run on B60's 16 identified coordinates, ON TOP of
B59b's route (kinetic_core_b59b_prereg.md Amendment 2); B59a-r on the old trunk is deferred until B59b's verdict.
