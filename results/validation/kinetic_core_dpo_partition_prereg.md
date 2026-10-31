# Pre-registration (DRAFT): the deoxyosone partition and the cysteine fate (plan step B7)

*Drafted 2026-10-31 ~02:00 from a read-only groundwork pass (no engine runs; the residuals below are read from the
unconverged fit records `s042_e` / `s042_r3en`, so they are diagnosis, not verdicts -- STATUS_AND_PLAN sec. 5 rule (a)).
Not yet reviewed; no fit runs before an independent design review. Experiment register EX-031.*

## 1. Why this outranks further MFT-route arms

The thiols are made downstream of two pools whose fate the engine gets wrong in the fit's own data: the
1-deoxypentosone (DPO) and 3-deoxypentosone (TDP) split, and free cysteine. Every MFT / FFT route arm so far inherits
both errors (thiol prereg secs. 15 - 16).

## 2. What the groundwork found

**The measurement model hides the partition error.** `series_of` keys calibration scales by paper x figure x species,
so TDP and DPO of one figure get separate scales. Those scales absorb 85 - 93 % of the TDP / DPO ratio error: 1.57 of
1.85 dex in Zhai 2019 Fig. 3b, 1.45 of 1.56 in Zhai 2023 Food Chem Fig. 2, 1.31 of 1.40 in Zhai 2023 JAFC Fig. 5bc,
0.66 of 0.72 in Zhai 2021 Fig. 5ab. The tau 0.3 shrinkage does not hold them: with 12 - 47 rows and fitted series noise
of 0.40 - 0.57 dex, each offset keeps 0.84 - 0.95 of the series' mean residual. Both species are quinoxalines read by
one method in one figure, so a shared scale with a response-factor prior is the defensible model.

**TDP has no real exit; DPO is over-drained.** TDP's absorbed offset is positive in every Zhai figure and in the
amine-free xylose pot (Zhou 2025b, +1.43); TDP lives ~35 h in the engine at 100 C (TDP -> furfural ~2e-6 /min, the
shared osone decay 4.8e-4 /min) against DPO's ~6 min. The fitted C2 + C3 drain on DPO (+2.42 dex) is the fit making
the measured absolute MGO out of DPO, the only MGO parent switched on (`r_tdp_mgo`, `r_pent_mgo` are off in
`trial()`) -- a compensation, and the DPO scale then hides it.

**The split's basis is one paper at one temperature.** CR-12's `dpo_share` 0.30 rests on Zhai 2019's fed-TTCA ratio
(120 C only), and a concentration ratio is a formation split only if the two lifetimes are equal (they differ ~300x
in the engine). The same laboratory measures TDP / DPO 17 - 28 at pH 5.5, 14 - 22 at pH 6, 6 - 8 at pH 7, 0.7 - 1.9 at
pH 8 (Zhai 2021, 100 C) and 6 - 15 at pH 7 (Zhai 2023): a strong pH dependence the engine's split does not have.

**The cysteine misfit is probably a ring-opening rate, not a missing sink.** Free cysteine is high (Zhai 2021 Fig. 6
+1.55 / +1.44 / +1.05 / +0.72 dex at pH 5.5 / 6 / 7 / 8; Zhai 2023 FC +0.87 / +0.33) while TTCA is LOW in the same pots
(Zhai 2021 Fig. 4b -0.59 at pH 5.5 to -0.14 at pH 8; Fig. 3a -0.83 / -1.41 at 120 / 140 C); the Cao and Zhao blanks
fit (-0.03 / -0.08). The engine's ring-opening (~0.028 /min, fitted, frozen since B2.3, no pH factor) is ~10x faster
than Zhai's total TTCA loss (0.0027 /min at pH 5.5); Cao's 0.4 M pot relaxes ~9x faster than k_open, which is why it
equilibrates within 30 min whatever k_open is. A slower ring-opening reconciles Cao with Zhai and Kang; TTCA-low with
Cys-high is that rival's signature (a pure extra sink would push TTCA lower still). Kang's 16 % ceiling still needs a
cysteine loss in TTCA pots ~10x faster than in Kang's cysteine-alone pot. **No Cys, TTCA or Amadori row is in the
joint cost today** -- the route table is scored only in dry checks.

## 3. Scale-free signatures (each within one figure; declared now, before any B7 arm)

| # | signature | measured | engine at R3En (unconverged) | tolerance |
|---|---|---|---|---|
| i | time slope of log10(TDP / DPO), Zhai 2019 Fig. 3b | +0.11 dex per decade of time | ~+1.5 (20 -> 297) | +-0.3 |
| ii | pH slope of log10(TDP / DPO), Zhai 2021 (both ladders) | -0.49 dex per pH unit | ~0 (no pH in the split) | +-0.15 |
| iii | cysteine double ratio, Zhai 2019: (TDP / DPO) with Cys over without | 0.46 (0.31 - 0.60) | to be read | [0.3, 0.7] |
| iv | free Cys pH 8 / pH 5.5 at 120 min, Zhai 2021 Fig. 6 | x4.8 | to be read | +-0.2 dex |
| v | Nedvidek 1992 glucose pH 5 -> 7: fragments / deoxyosones 2:3 -> 7:1; 1-DG / 3-DG x5 | direction | to be read | sign |

## 4. Rivals (all under ONE new measurement model)

**Measurement-model change, applied to every arm including the null:** (1) TDP and DPO of one figure share one scale,
with a declared response-factor prior N(0, 0.2^2) dex on their difference; (2) the route table's Cys, TTCA and
Amadori-compound rows enter the cost. The null therefore tests the CURRENT structure under the honest model.

| arm | structure | added coordinates |
|---|---|---|
| N | R3En / R4f structure (whichever the converged comparison admits) | 0 |
| R-A | pH-dependent enolisation split (de Roos 2005 acid / base; Nedvidek 1992; Zhai 2021) + TTCA ring-opening as a coordinate with a base slope | `dpo_share` at pH 7, split pH slope, `offset_k_ttca_cys`, ring-opening pH slope (4) |
| R-B | cysteine + aldehydic dicarbonyl (TDP, glyoxal, MGO) -> thiazolidine, one class constant with a free-amine factor (the family of `r_cys_actz`, `r_hmf_cys`; Zhai 2019 Fig. 3b's cysteine suppression of TDP 0.20, GO 0.23 vs DPO 0.45) + TDP's own decay split from `k_osone_decay` (3-DX -> formic + C4, by analogy with Martins' measured 3-DG -> formic) | 2 |
| R-AB | R-A + R-B | 6 |

R-B's literature support is weak (searches found patents and a penicillamine analogue, no rate constant): if the
data choose it, the constant is a fitted one with that label. The MGO parent question (c2c3) is not addressed by
any arm; it is the diagnostic below.

## 5. Decision order (declared now)

1. Converged fits only (rule a). 2. Signature (i). 3. Signature (ii). 4. Signature (iii). 5. Zhai 2021 Cys AND TTCA
both within 0.3 dex, and signature (iv). 6. Must-not-break: Cao and Zhao blanks within 0.1 dex, Kang <= 16.3 % free
cysteine, Zhai MGO absolute, Cerny pot D (no worse than the null). 7. Cost under the fixed measurement model.
8. Diagnostic: `offset_k_dpo_c2c3` below +1.5.

**Expected (written now, before any arm):** R-A passes (ii) (P 0.6); R-A brings Zhai 2021's Cys and TTCA within
0.3 dex (P 0.5); R-B passes (iii) (P 0.45); c2c3 drops below +1.5 in any arm (P 0.3 -- MGO still needs a parent);
the null fails (i) and (ii) under the new measurement model (P 0.85).

## 6. Open before the design is final

- The shared-scale change alters every earlier arm's cost: the round-3 / R4f comparisons stay under the old model
  (one objective each); B7's arms are compared only with each other.
- Implementation touches `series_of` / `method_of` and the fit tables in the joint fit script, and the CAMA split
  in `parameters_sulfur.py` -- not before the fit queue finishes (a running fit imports both).
- An independent design review, then a premortem of B7 against the 48 h time-box of section B.

## 7. Independent design review (2026-10-31 ~02:20, no engine runs) and the revised scope

**Checked against the first converged fit (s042_e_conv):** the scale absorption holds (85 / 90 / 94 / 93 % of the
TDP / DPO ratio error in Zhai 2019 / 2021 / 2023 FC / 2023 JAFC; magnitudes 0.07 - 0.27 dex below sec. 2's, which were
R3En's); TDP high in every figure holds (+0.40 to +1.03; Zhou 2025b +1.43); "DPO over-drained" is arm-dependent (in
e_conv DPO is roughly right because its DDP route sits at the lower bound); c2c3 +2.08. Cys-high / TTCA-low cannot be
checked from the fit records (the route rows are in none).

**BLOCKING.** (B1) The route rows would dominate the cost: 357 rows (96 in Zhai 2021 Fig. 6), printed SDs median 0.015
dex, absolute, no level term -- the cysteine bias alone would cost ~2,000 - 6,000 units against e_conv's 507, turning
B7 into a cysteine fit; and --noise-from would refuse (no noise for those series). (B2) A slower ring-opening cannot
reconcile Cao with Zhai: Cao's free Cys falls 400 -> 75 mM by 30 min (relaxation >= 0.12 /min, so k_open(90 C) >=
0.013 /min), while Zhai 2021 at 100 C caps k_open <= 0.0033 /min (nearer 4e-4): a 4 - 30x contradiction, worse with any
Ea -- sec. 2's "whatever k_open is" is wrong. Cao and Zhao use 0.2 M phosphate, Zhai water: a missing rival R-P,
general-acid (buffer) catalysis of the ring-chain step and of r_cama_form (de Roos 2005 calls the thiazolidine ->
Amadori step general-acid catalysed; r_cama_form's anchor is Cao's phosphate pot). (B3) The TDP exit belongs in its
own arm: Zhai 2023 JAFC's TDP falls after its peak at an apparent 0.009 /min (120 C) and 0.025 /min (140 C), against
the engine's ~5e-4.

**SHOULD-FIX.** (S1) The shared scale needs code: MeasurementModel has one level per row -- use level_of for the
figure scale, keep series_of per species (pooled_noise per species), add the response-factor term in fitting.py with
tests. Same method within each figure is confirmed (OPD quinoxalines, HPLC-DAD), but calibration differs by paper and
no authentic 3-DX / 1-DX standard is reported, so the 0.2 dex prior is unsourced: sensitivity at 0.4; the shared scale
applies to the Zhou and Desclaux figures too. (S2) Signatures recomputed and confirmed ((i) +0.09 to +0.11 dex /
decade; (ii) -0.51 / -0.46 dex / pH; (iii) 0.45; (iv) x4.8, x5.3), but (i) is scored per figure (other figures -0.52 to
+1.11), (ii) is curved (segment slopes -0.14 / -0.21, -0.38, -0.70 / -0.80) and is scored by segment, (iii) is
contradicted by Zhai 2023 JAFC (0.53 - 2.6: a conflict to register), (v) cannot be computed in the sulfur lane (prior
only); add (vi) TDP's post-peak decline; all are in-sample, with a chance baseline. (S3) "Not a missing sink" is
premature: at pH 5.5 / 6 TDP + DPO (0.65 / 0.99 mM) exceed free Cys (0.33 / 0.43 mM), so some cysteine stays bound or
is lost; the null takes R4f's k_ttca_deg. (S4) Split R-A into the enolisation-split law and the ring-opening law; with
detailed balance at fixed K a ring-opening pH law changes the rate, not the equilibrium. (S5) GO is not a sulfur-lane
species: drop it from R-B or add it. NOTES: no rule if every arm fails (i); the v2 MFT groups join the must-not-break
list; P(R-A breaks Cao) ~0.7; the "engine" column is re-read at the null; the converged R3En is not on disk yet.

**Revised scope (adopted).**
- **B7a (inside the 48 h box): the deoxyosone partition only, deoxyosone table only.** Measurement model: figure level
  (level_of) + per-species series with a response-factor prior (0.2 dex; sensitivity 0.4), coded in fitting.py with
  tests, noise frozen once under it. Arms: N (the admitted structure + R4f's k_ttca_deg); R-S (an enolisation-split pH
  law: dpo_share at pH 7 + a slope or a two-segment law, 2 coordinates); R-T (a TDP exit: 3-DX -> formic acid at
  Martins 2005b's k5 as a declared transfer, 0 coordinates, plus the coded base-catalysed r_tdp_mgo switched on -- a
  second MGO parent so the c2c3 diagnostic becomes testable). Signatures (i) per figure, (ii) by segment, (vi);
  (iii) reported, not decisive. ~6 h of work + ~3 h of compute.
- **B7b (only after a no-fit check of the null's cost by block): the cysteine fate** -- route rows with an
  between-paper level term and a noise floor (or as signatures / must-not-break, capped per series), rivals R-P
  (buffer catalysis of the ring-chain step and r_cama_form) and the ring-opening law; otherwise backlog.
- **Cut:** R-AB, signature (v) as a score.
- **Expected (written now, before any B7a arm):** R-S passes (ii) by segment (P 0.45); R-T passes (vi) (P 0.5) and
  lowers c2c3 below +1.5 (P 0.35); the null fails (i) in >= 2 figures (P 0.8).

## 8. B7a arms, specified (2026-10-31 ~03:50, after the null's signatures at R3En (EX-033), before any B7a arm)

**Null at the admitted structure (EX-033, one-level measurement model):** time slope 1 / 8 pots (the engine's
log10(TDP / DPO) climbs +1.1 to +1.8 per decade where the measured is flat), pH segments 1 / 6 (the engine's ratio
RISES with pH, the measured falls -0.14 to -0.80 per pH unit), TDP decline x8 too slow at 120 C (right at 140 C).
The B7a null under the two-level measurement model (`--within-tau` 0.2 and 0.4, R3En from `s042_r3en_fixed`, noise
re-estimated) is running; its signatures are read with the same scorer.

**All arms:** R3En's structure and coordinates (R4f showed `k_ttca_deg` stays near its value, so it is not freed),
the two-level measurement model at within_tau 0.2 (0.4 as the sensitivity of the winner), series noise taken from
the B7a null at the same width (one objective), from the B7a null's optimum, budget 1200, ftol 1e-4.

- **N** = the B7a null itself.
- **R-S (the enolisation split's pH law), +2 coordinates.** `CysAmadoriRoute.dpo_share` becomes
  logistic(logit(s7) + b (pH - 7)) with s7 the 1-deoxyosone share at pH 7 (kind acid_yield, [0.02, 0.98], start
  0.30 = CR-12) and b in logit units per in-situ pH unit ([-3, +3], start 0). Basis: de Roos 2005 (acid -> 1,2-
  enolisation / 3-deoxyosone, base -> 2,3-enolisation / 1-deoxyosone), Nedvidek 1992 (1-DG / 3-DG x5 from pH 5 to 7),
  Zhai 2021's ratio falling with pH. With b = 0 and s7 = 0.30 the engine is exactly today's (short-circuit, no
  logit round trip).
- **R-T (a TDP exit), 0 coordinates.** TDP's own decay (`r_osone_decay_tdp`, today the shared fitted osone decay
  4.8e-4 /min) at Martins 2005b's 3-DG -> formic acid step k5 as a DECLARED TRANSFER (hexose 3-deoxyosone ->
  pentose): 1.0e-2 /min at 100 C (geometric mean of the four non-zero cells 6.3e-3, 2.0e-2, 4.1e-3, 1.9e-2 /min,
  pH 5.5 - 7.5; Martins reports "low pH-dependence", Table 3 pD 0.09 +- 1.10, so no pH term), Ea 29.6 kJ/mol (Martins
  2005, step 5). Products unchanged (ACID + FRAG_C 4: the formic acid joins the lumped acid pool the charge balance
  already carries). Implemented as a derived rate behind a switch (off: bit-identical). Arithmetic check: 1.6e-2 /min
  at 120 C, 2.5e-2 at 140 C vs Zhai 2023 JAFC's apparent post-peak declines 0.0097 / 0.0266.
- **R-ST = R-S + R-T**, +2 coordinates (the null's two failures are in different signatures -- time and pH -- which
  the two arms address separately, so their combination is declared now rather than added after a partial result).

**Decision order (sec. 5, as revised in sec. 7):** converged fits only; (i) per pot (pass >= 5 / 8); (ii) by segment
(pass >= 4 / 6); (vi) (2 / 2); must-not-break: the v2 thiol groups no worse than the null, Cao / Zhao blanks, Zhai
MGO absolute, Cerny pot D no worse than the null; cost under the fixed measurement model (R-T has no new coordinate,
so its cost difference is a direct test; R-S / R-ST pay 2 coordinates); diagnostic: c2c3 below +1.5.
**Expected (written now):** R-T passes (vi) (P 0.6) and >= 5 / 8 of (i) (P 0.35); R-S passes >= 4 / 6 of (ii) (P 0.35
-- the measured pH law is curved and a logit-linear one may miss a segment); R-ST passes all three families (P 0.2);
R-T lowers the cost below N's (P 0.6); c2c3 below +1.5 in any arm (P 0.25 -- MGO still needs a parent).

## 9. Design review of sec. 8 (2026-10-31 ~04:00; no engine runs) and the revised arms

**BLOCKING.** (B1) R-T's transfer was a forking-path choice: the Ea 29.6 is Martins 2005's, whose own constant for
the step is 3.5e-2 /min -- the trunk's `k_tdg_fa` 3.45e-2 -- while 2005b's cells scatter ~0.35 dex (its pH 6.8 cell
8.5x lower) and the Ea itself is contested (Knol 84 +- 14, Yu 43.9; at 84 the rate fails (vi) by 0.6 - 0.7 dex); the
source was picked after the declines were read. (B2) At 140 C R-T SLOWS TDP's sink: the null's shared osone decay (B60,
Ea 143.5) gives 4.7e-4 / 5.0e-3 / 4.2e-2 /min at 100 / 120 / 140 C, so R-T's 0.025 replaces 0.042 exactly where the
null passes both signatures; an apparent post-peak decline is only a lower bound on a first-order exit while
formation continues. (B3) A logit-linear split gives a constant log-ratio slope against measured segments -0.21 /
-0.38 / -0.80 and -0.14 / -0.38 / -0.70: at most 2 per ladder can pass, so R-S's fail was close to predetermined;
and sec. 7's base-catalysed `r_tdp_mgo` (a sink with the right pH direction and a second MGO parent) was silently
dropped. SHOULD-FIX: (S1) the declines were selected -- the JAFC 100 C pot (~0.018 /min) and the +Cys pots (0.025 -
0.035, flat in temperature) were left out, Zhai 2023 FC's 140 C pot plateaus, the apparent declines do not rise with
temperature (no single Arrhenius exit), and **in JAFC TDP and DPO fall at the same rate (140 C: 0.025 vs 0.028 /min)
with the ratio flat from 10 min -- equal lifetimes**: the missing rival is a LONG-LIVED DPO (the c2c3 drain released)
plus an MGO parent. (S2) No chance baseline for these signatures, and the thresholds were set after reading the null.
(S3) The split acts only on the Amadori flux (r_pent_caramel, r_pent_thermal_dpo have no pH factor); the Zhai 2021 pots
are unbuffered, so the in-situ pH span is narrower than the labels. (S4) s7 is identified mainly by the response-
factor term; report cost by block and profile the law with s7 fixed; TDP / DPO are 276 of 521 rows and noise
re-estimated at N absorbs N's misfit (compare with e_conv's 0.40 - 0.57). (S5) R-T's extra ACID lowers the in-situ pH
of unbuffered pots: report the shift. Must-not-break additions: the null's two passes at 140 C, furfural (17
calibrated rows; r_tdp_fur scales with TDP) and FFT, Zhou 2025b's amine-free TDP, Kang <= 16.3 %, the hold-outs.

**Revised arms (adopted; all else as sec. 8).**
- **B7S** -- a CURVED split law with the same two coordinates: TDP : DPO odds(pH) = odds7 (1 + 10^(7 - pHc)) /
  (1 + 10^(pH - pHc)) (base-saturating 2,3-enolisation; s7 = DPO share at pH 7, pHc in [5, 10]).
- **B7T** -- TDP's own decay at the trunk's own hexose pairing (`k_tdg_fa`, its declared k_ref / t_ref / Ea),
  band +-1 dex declared, sensitivities at Martins 2005b's 1.0e-2 /min and at Ea 84 (Knol); 0 coordinates.
- **B7ST** = B7S + B7T.
- **B7D (new, the review's rival)** -- a long-lived DPO with a second MGO parent: the base-catalysed `r_tdp_mgo`
  switched on with one coordinate (log10 k at 100 C / pH 5, [-8, -1], start -5.08 = sub-model 3 arm A3's fit; Ea at
  its declared 100 kJ/mol), `offset_k_dpo_c2c3` free as before (the data can lower it once MGO has another parent).
- Scorer: all post-peak declines reported (the two JAFC TTCA pots decisive), DPO's decline beside TDP's, Zhai 2021's
  in-situ pH; a chance baseline for these signatures before any arm is credited.
**Expected (rewritten, post-review):** B7T passes (vi) 2 / 2 (P 0.35) and >= 5 / 8 of (i) (P <= 0.15) and loses at
least one of the null's 140 C passes (P 0.5); B7S passes >= 4 / 6 segments (P 0.35); **B7D passes >= 5 / 8 of (i)
(P 0.35) and brings c2c3 below +1.5 (P 0.4)**; B7ST all three families (P 0.15).

## 10. Before any B7a arm is read (2026-10-31 ~04:15): the implementation, its review, and two declared changes

Implemented on branch `b7a-two-level-measurement` (B7S / B7T / B7ST / B7D; every default bit-identical -- 1048 rate
constants equal to the parent commit with the switches off; 154 targeted tests; the switch, choices, unpinned-constant
and evidence-anchor gates pass) and code-reviewed independently (nothing blocking). **Two starts changed from sec. 9,
declared now, before any arm is fitted or read:** B7S starts at pH_c 9.5 (within ~3 % of the null at pH 8; at 7.5 the
start already moved the split, so the arm would not have started at the null); B7D starts at log10 k_tdp_mgo = -7
(sub-model 3 arm A3's -5.08 was fitted with its barrier at the 40 kJ/mol bound; at the declared 100 kJ/mol it gives
~1500x the null's MGO at 120 - 140 C, far from the null). B7T's primary constant is the trunk's `k_tdg_fa` (3.45e-2 /min
at 100 C, Ea 30), ~5.6e-2 at 120 C against Zhai's apparent 0.0097: its (vi) fail at 120 C is likely (P 0.6) and is
part of the test, not a reason to swap the constant (the Martins 2005b 1.0e-2 sensitivity runs beside it). A chance
baseline of these signatures (`score_deoxyosone_signatures.py --chance`) runs before any arm is credited.

## 11. B7a arms: results (2026-10-31 06:20) -- no arm passes; the diagnosis moves to DPO's pH-dependent sinks

All four converged (status 2) under the B7a null's fixed series noise (w 0.2), from the null's optimum.

| arm | cost | added | key constants | (i) time | (ii) pH segments | (vi) decline | thiol groups | Cerny pot D |
|---|---|---|---|---|---|---|---|---|
| B7D | 1074.4 | 1 | log10 k_tdp_mgo -6.75 +- 0.19; c2c3 +2.37 | 2 / 8 | 0 / 6 | 1 / 2 | 3 / 11 | 0.992 |
| B7S | 1038.9 | 2 | DPO share at pH 7 0.62 +- 0.13; pH_c 6.79 +- 0.41; c2c3 +2.32 | 1 / 8 | 0 / 6 | 1 / 2 | 2 / 11 | 0.994 |
| B7T | 1146.0 | 0 | (trunk k_tdg_fa) c2c3 +2.45 | 1 / 8 | 1 / 6 | 2 / 2 | 4 / 11 | 0.993 |
| B7ST | 1082.0 | 2 | share 0.69 +- 0.13; pH_c 6.50 +- 0.36 | 2 / 8 | 0 / 6 | 2 / 2 | 3 / 11 | 0.996 |

The null's own cost under the same fixed noise (`s042_r3en_b7anull_w02_fixed`, converged) is **1093**: B7S -54 for two
coordinates and B7D -19 for one are real cost gains (delta chi-square / reduced chi-square ~24 and ~8); B7T +53 for none is
rejected by the fit rows; B7ST -11 for two is weak. **Against the declared thresholds -- (i) >= 5 / 8, (ii)
>= 4 / 6, (vi) 2 / 2 -- no arm passes.** Expectations: B7D (i) (P 0.35) no, c2c3 < 1.5 (P 0.4) no; B7S (ii) (P 0.35)
no; B7T (vi) 2 / 2 (P 0.35) YES, (i) (P <= 0.15) no; B7ST all three (P 0.15) no. B7T's 3-DX exit at the trunk's
hexose constant reproduces the post-peak declines but RAISES the cost by ~44 for no coordinate: the fit rows reject the
exit at that rate. B7S lowers the cost by ~63 for two coordinates, but its pH law does not show in the ratio.

**Why the split's pH law does not reach the ratio (flux attribution at B7S's optimum, no fit).** DPO's sinks are
pH-dependent the OTHER way: at pH 8 the base-favoured beta-cleavage to acetic acid (`r_dpo_aa`, slope fitted on
Davidek 2008's glycine pots) takes 51 - 61 % of DPO, at pH 5.5 the DDP route takes 58 %; so DPO is drained hardest
exactly where the split makes more of it. And 39 - 45 % of TDP comes from the amine-free pentose caramelisation
(`r_pent_caramel`), which has no pH law. DPO is made mostly through the Amadori compound (56 - 92 %), so the split
does act on it -- the sinks cancel it.

**What this points at (for the owner's review; no arm is pre-registered yet).** A tension between Davidek 2008
(xylose + glycine: acetic acid from DPO strongly base-favoured) and Zhai 2021 (TTCA, no free amine but cysteine: DPO
relatively MORE abundant at high pH). Candidate: DPO's beta-cleavage to acetic acid is AMINE-ASSISTED (glycine in
Davidek's pots, none in Zhai's), so its base factor should multiply a free-amine term, not act alone -- one structural
change, no new coordinate if the amine order is declared. Second: TDP's caramelisation source and its pH dependence
(Zhai 2021's xylose ladder isolates it). Both are testable on data already keyed; neither was run tonight.

## 12. Notes for the cysteine fate (B7b, folded into the hub step) from the new dossiers (2026-10-31 ~08:30)

Roberts 1987 (roberts1987_extraction.md): no ring-opening rate, K or buffer data; the sugar-cysteine thiazolidines
anomerise within minutes in unbuffered water at room temperature -- mild evidence AGAINST an intrinsically slow ring
opening in water, so buffer catalysis is not supported as the whole explanation of Cao (phosphate) vs Zhai (water);
leads for a rate: Nagasawa 1981 (10.1002/jhet.5570180538), Schubert 1939 (10.1016/s0021-9258(18)73532-1). Parker 2014:
one 140 C / 30 min end point -- MFT goes the opposite way to norfuraneol across ribose vs ascorbic acid (consistent
with Cerny's minor norfuraneol route); nothing on the 100 C lag. Tai & Ho 1997 / 1998: no free thiol, H2S or MFT;
sulfur volatiles carry < 0.5 % of the charged sulfur, so the missing cysteine sink is non-volatile. **A rival worth
testing:** headspace O2 oxidises thiol to disulfide in proportion to the headspace, not the charge -- large in the
5 - 10 mM Zhai pots, <= 10 % at Cao's 0.4 M, which matches where the engine misfits; against it, Zhai 2021's cysteine
excess is largest at pH 5.5 (thiolate oxidation predicts the opposite). P(oxidation is the main cause) ~0.2; the Zhai
vessels' headspace is the deciding fact (check the dossiers).

## 13. B7a-2, revised by the Davidek 2006a / 2006b / Novotny 2008 dossiers -- design, 2026-10-31 ~09:30, before any fit

**What the dossiers changed.** The amine-assisted cleavage (sec. 11's candidate) is contradicted: synthetic
1-deoxyglucosone heated WITHOUT amine in 0.16 - 0.2 M phosphate gives 0.80 - 0.84 mol acetic acid / mol at pH 6.5 - 8
(60 min, 120 C; Davidek 2006a Fig. 3), 0.56 at pH 6, 0.20 at pH 5 -- an amine REQUIREMENT is excluded for the hexose
analogue (P 0.85; 0.75 that it transfers to the pentose), amine ACCELERATION is untested (no paired run). The pH law
is a SATURATING step, flat above pH 6.5 ("no significant differences ... at pH 6.5, 7, 8"), not the engine's
log-linear base factor 10^(aa_ph_slope (pH - 6)), which keeps raising DPO's beta-cleavage from pH 7 to 8 -- exactly the
over-drain of DPO at high pH that sec. 11 found cancelling the split law. Phosphate cannot be the Davidek / Zhai
difference (Zhai 2021's TTCA pots are also 0.2 M phosphate). The 1-DG yields are branching fractions (1-DG is gone
within 5 min at pH >= 6), so the curve constrains the cleavage's share of DPO's fate, not its rate directly.

**Arms (both under the B7a null's fixed series noise, w 0.2, from the null's optimum; budget 1200, ftol 1e-4):**
- **B7AS -- the saturating law alone.** r_dpo_aa = k_max / (1 + 10^(pK - pH_in_situ)) (an acid-base step completing
  above pK; the rate saturates), replacing 10^(aa_ph_slope (pH - 6)). Coordinates: log10 k_max (replacing log10
  k_dpo_aa at pH 6) and pK in [5.0, 7.5] (replacing aa_ph_slope) -- the SAME number of coordinates as the null, so
  the cost difference is a direct test. Reference band from Davidek 2006a with a pH-flat competing sink: pK ~5.8 -
  6.4 (reported against, not imposed).
- **B7ASS -- the saturating law + B7S's split law** (+2 coordinates, the split's share at pH 7 and pH_c).
**Decided by** (sec. 7's order): converged; the partition's scale-free signatures against their chance baseline
(EX-037: time slope >= 5 / 8 at chance 0.02; pH segments >= 4 / 6 at chance 0.00); must-not-break: the Davidek 2008
acetic plateaus and pH 8 / 6 ratio (the rows the slope was fitted on), the v2 thiol groups no worse than the null;
cost vs the null under the same noise (1093) and vs B7S (1039).
**Expected (written now):** B7AS's fitted pK within 5.8 - 6.4 (P 0.4); B7AS passes >= 4 / 6 pH segments (P 0.15);
B7ASS passes >= 4 / 6 (P 0.3) and >= 5 / 8 time slopes (P 0.1 -- the time slope is a lifetime problem the pH law
does not touch); the Davidek 2008 acetic signatures hold in both (P 0.6). If neither arm moves the pH segments,
section B closes by its rule (24 h cap, 11-01 ~08:00).
**Note at implementation (before any fit):** with pK 6 the saturating factor is 0.18 / 1 / 1.52 / 1.82 / 1.98 at pH 5 /
6 / 6.5 / 7 / 8, while R3En's fitted slope (0.15 dex per pH unit) gives 0.71 / 1 / 1.19 / 1.41 / 2.0: the two laws
nearly coincide from pH 6 to 8 and differ mainly BELOW pH 6. So B7AS mostly tests the acid side; it cannot by itself
remove a high-pH over-drain that the slope did not cause. The P 0.15 for B7AS's pH segments stands.
**Review of B7a-2 (independent, before any fit; nothing blocking) -- amendments adopted now.** (1) The coordinate is
log10 k at pH 6 (the saturating factor is normalised to 1 at pH 6), not log10 k_max as written above. (2) Davidek's
5.8 - 6.4 band is on the LABEL scale (set at room temperature); the engine's pK acts on the in-situ pH, ~0.2 - 0.3 lower
for phosphate at 120 C: the comparison band is ~5.5 - 6.2. (3) The pK test is not independent: the pentose fit rows are
pH 6 and 8 only, both laws are one-parameter maps of the same 8 / 6 ratio, and slope 0.15 maps to pK ~6.0 -- a pK in
band is not evidence. (4) The informative check is the BRANCHING share: the engine's DPO -> acetic acid flux share at
pH 5 / 6 / 6.5 / 7 / 8 against Fig. 3's shape (flat above 6.5) -- reported for both arms. (5) Cost: with the same
coordinate count, B7AS counts as a gain only if it is >= 5 units below the null under the same noise (delta chi-square
>= 10, ~2 after the reduced chi-square); B7ASS against B7S (1039) with the same rule. (6) The chance baselines (EX-037)
were drawn in R3En's box; replacing one coordinate by another of the same reach is taken as negligible (stated, not
re-drawn). **Expectations revised:** B7AS >= 4 / 6 pH segments P 0.07; B7ASS P 0.1; the Davidek 2008 acetic
signatures hold P 0.75; the engine's DPO -> acetic share flat (within 0.1) from pH 6.5 to 8 under B7AS P 0.5.

**CORRECTION (2026-10-31 ~08:45, before the B7a-2 arms were read).** Sec. 13 states "Phosphate cannot be the Davidek /
Zhai difference (Zhai 2021's TTCA pots are also 0.2 M phosphate)". That is WRONG: Zhai 2021's TTCA pots are
UNBUFFERED water, pH set with NaOH / HCl (zhai2021_extraction.md, Methods; the fit table's pots), as are Zhai 2019,
Zhai 2023 and Kang 2026. The claim came from a dossier agent's summary and was copied without checking. Consequences:
Davidek 2006a's amine-free 1-DG pots (0.16 - 0.2 M phosphate) and Zhai's TTCA pots (water) differ in phosphate, so
**phosphate (general-acid / base) catalysis of the DPO beta-cleavage is a live rival** for the Davidek / Zhai
tension -- the cysteine mechanism review (scratchpad mechanism_reviews/cysteine.md) independently ranks catalysis by
solutes and buffer as the leading cause of the cysteine misfit (P 0.55), with the +-xylose contrast of Zhai 2021
Fig. 6 (free cysteine RISES on adding xylose: no equilibrium model allows it) as its on-disk discriminator. The B7a-2
arms (saturating pH law) are unaffected as tests; their interpretation will name phosphate as the untested rival.

## 14. A held-out signature declared before the B7a-2 arms are read (2026-10-31 ~08:50), and the hub review's leads

**Signature (vii), declared now (held out from every B7 arm's credit until a fit is designed on it):** the +xylose /
-xylose double ratio of TDP / DPO by pH in Zhai 2021 Fig. 5 (100 C, 120 min), measured 0.85 / 0.92 / 1.22 / 2.7 at pH
5.5 / 6 / 7 / 8 (the hub mechanism review's arithmetic on the keyed tables): extra xylose re-traps released cysteine,
and 3-DX is rescued only at high pH -- the signature of a cysteine trap of the aldehydic 3-DX that grows with pH.
Reported for every arm from B7a-2 on; tolerance +-0.15 dex per point, pass = 3 of 4.

**The hub review's ranked hypotheses (scratchpad mechanism_reviews/hub.md; inputs to the hub step, roadmap step 5):**
H1 both deoxyosones quasi-steady (short-lived vs the cysteine Amadori compound) while the engine's 3-DX lives ~35 h
(P 0.65; Hofmann 1999, Zhai 2023 JAFC, Zhai 2021 pH 5.5, Martins M4); H2 free cysteine traps 3-DX faster at higher pH
(P 0.45; signature vii); H3 the split itself moves toward 1-DX with pH (P 0.5 as a contributor; Martins 2005b pD,
Nedvidek 1992); H4 the 1-DX sinks calibrated in phosphate and applied to water (P 0.35; see the correction above);
H5 the pH-flat caramelisation source of 3-DX too large (P 0.4); H6 O2 diverting the Amadori compound to the osone (P
0.15; Gobert 2009). Null for the hub step: Martins & van Boekel 2005 M4 (its topology already puts both osones in
quasi-steady state).
(Timing note: the command that appended this section printed the B7a-2 log first, so B7AS's COST -- 1079 -> 1022 --
was seen just before (vii) was written; none of B7AS's signature values had been computed or read.)

## 15. B7a-2 results (2026-10-31 09:20) and the close of section B's partition work

| arm | cost (null 1093) | added | key constants | (i) time | (ii) pH segments | (vii) +xyl | Cerny pot D |
|---|---|---|---|---|---|---|---|
| B7AS | 1021.6 | 0 | pK 5.35 +- 0.32 | 1 / 8 | 0 / 6 | 2 / 4 | 0.994 |
| B7ASS | 901.5 | 2 | pK 5.64; DPO share at pH 7 0.87; pH_c 6.38 | 0 / 8 | 1 / 6 | **3 / 4 (pass)** | 0.997 |

Both converged (status 3). **Neither passes (i) or (ii); B7ASS is the first arm to pass (vii).** Under the owner's rule
(mechanisms before constants) the large cost gains (-71, -192) are not evidence. What B7ASS does show: the split law
together with the saturating beta-cleavage turns the ratio's pH slope to the RIGHT SIGN from pH 6 to 8 (-0.23 / -0.15,
-0.22 / -0.10 dex per pH unit vs measured -0.38 / -0.80, -0.38 / -0.70; too shallow) and reproduces the +xylose
double ratio at pH 5.5 - 7 -- while the time slope gets WORSE (+0.4 to +2.0 per decade vs ~flat). So the pH structure
and the time structure are separate faults: the time slope is the 3-DX LIFETIME (the hub review's H1, P 0.65: both
deoxyosones quasi-steady, the engine's 3-DX ~35 h), which no arm here touched. Expectations: B7AS pH segments (P 0.07)
no; B7ASS (P 0.1) no; the B7AS pK in the in-situ band ~5.5 - 6.2 -- 5.35 no, B7ASS's 5.64 yes.

**Section B's partition work closes here** (the 24 h cap, roadmap row 1). Carried to the hub step (roadmap row 5),
each to be tested for FEASIBILITY before any fit (STATUS_AND_PLAN, mechanisms before constants): (a) 3-DX
quasi-steady (Martins M4's lifetimes as the null); (b) the split law + saturating beta-cleavage (B7ASS's pH structure);
(c) a cysteine trap of 3-DX growing with pH (signature vii at pH 8); (d) catalysed TTCA hydrolysis (the cysteine
review); (e) phosphate as a variable (Davidek in phosphate, Zhai in water). Decisive tests: (i), (ii), (vii), the
+-xylose free-cysteine contrast, Cerny pot D.

## 16. The first hub FEASIBILITY test (mechanisms before constants) -- design, 2026-10-31 ~09:40, before any run

**Question.** Is there ANY constant set, inside declared bands, for which the structure "B7ASS's pH structure + a
real 3-DX exit" (arm B7ASST: the split's base-saturating law, the saturating beta-cleavage, and TDP's own decay)
reproduces the partition's scale-free signatures? No fit: a structure that cannot reach them anywhere in its bands is
rejected; one that can goes to a pre-registered fit.

**Sampled (uniform, independent), all else held at B7ASS's optimum (s042_b7ass_w02):** dpo_share7 [0.02, 0.98];
split_ph_c [5, 10]; aa_pk [5, 7.5]; log10_k_dpo_aa (pH 6) its optimum +-1 dex; offset_k_dpo_c2c3 [0, 3] (the C2 + C3
drain, compensation-suspect); offset_k_dpo_ddp its optimum +-1; TDP's own decay log10 k(100 C) in [-3.5, -1] /min
(lifetimes ~5 min - 50 h; Martins' 3-DG ~30 min sits inside) and its Ea in [30, 100] kJ/mol. 300 draws (seed
20261031). **Null comparison:** the R3En-box chance baseline (EX-037) never reached (ii) >= 4 / 6 in 192 draws.
**Feasible** if any draw passes (i) >= 5 / 8 AND (ii) >= 4 / 6 AND (vii) >= 3 / 4; reported: the share of draws
passing each, the best draw per family, the joint passes, and the constants of the joint passes (are they physical?).
**Expected (written now):** feasible (>= 1 joint pass) P 0.3; (i) >= 5 / 8 reached by some draw P 0.55 (a 3-DX exit is
exactly what (i) needs); (ii) >= 4 / 6 reached P 0.35.

## 17. Feasibility v1 result, and v2 (the sources) -- 2026-10-31 ~10:05; v2 designed before it runs

**v1 (B7ASST, 299 of 300 draws; results/validation/kinetic_core_deoxyosone_feasibility_b7asst.md): NOT feasible.**
Time slope >= 5 / 8 reached by 0 draws (best 4); pH segments >= 4 / 6 by 6 (the R3En box: 0 of 192 -- this structure
CAN produce the pH shape); (vii) by 91; pH + (vii) together by 6; nothing with the time slope. Expectations: feasible
(P 0.3) no; (i) reachable (P 0.55) no; (ii) reachable (P 0.35) yes. **A 3-DX exit as fast as 0.1 /min cannot flatten
the ratio's time course**, so the time slope is not (only) 3-DX's lifetime: the SOURCES must shift in time -- the
hub review's H5 (pH-flat pentose caramelisation feeding 3-DX as free xylose is released from TTCA) and the early
1-DX burst from TTCA's direct breakdown fading (sec. 15's flux attribution).

**v2 (arm B7ASST2 = B7ASST + offset_k_pent_caramel + offset_k_ttca_deg), same tester, same reference and seed, 300
draws.** Added bands: offset_k_pent_caramel [-3, +0.5]; offset_k_ttca_deg [-3, +0.5]; offset_k_dpo_nf its optimum
+-1 (the tester's band table gained these three rows at v2; v1's record lists the bands it used). Feasible as in sec. 16.
**Expected (written now):** feasible P 0.25; time slope >= 5 / 8 reachable P 0.45; the joint passes, if any, need
caramelisation lowered by >= 1 dex (P 0.6 given feasible).

## 18. Feasibility v2 result, and v3 -- 2026-10-31 ~10:00; v3 designed before it runs

**v2 (B7ASST2, 299 draws; ..._feasibility_b7asst2.md): NOT feasible, but the time slope is now reachable.** Time slope
>= 5 / 8 by 3 draws (v1: 0), pH segments >= 4 / 6 by 4, (vii) by 64; time + pH together by 0, pH + (vii) by 2.
Expectations: feasible (P 0.25) no; (i) reachable (P 0.45) yes. **The sources drive the ratio's time course** (v1's
3-DX exit could not). Whether time and pH are in TENSION or the joint region is merely small in nine dimensions is
not decidable from 300 uniform draws.
**v3:** the same arm and bands, 600 draws, another seed (20261101), and the tester now records the sampled constants of
every draw reaching each family (the first 40) -- so the time-reaching and pH-reaching regions can be compared.
**Expected (written now):** a joint time + pH pass in 600 draws P 0.3; the time-reaching draws have caramelisation
lowered (offset < -1) in >= 2 of 3 cases P 0.6; if the two regions are disjoint in a constant (e.g. the split law's
pH_c), that constant's mechanism is the tension -- reported, then designed.

## 19. Feasibility v3: the tension located; v4 designed (2026-10-31 ~10:50, before v4 runs)

**v3 (B7ASST2, 597 of 600 draws, seed 20261101): NOT feasible, and the tension is located.** Time slope >= 5 / 8 by 3
draws, pH segments >= 4 / 6 by 11, (vii) by 110; time + pH jointly 0. The time-reaching and pH-reaching draws are
DISJOINT in the 3-DX exit: time-reaching log10 k(100 C) -1.13 to -1.01 (lifetime ~13 min, Ea ~91), pH-reaching -3.40
to -1.34 (slow). They also differ in the split (DPO share at pH 7: time 0.07 - 0.43, pH 0.34 - 0.94) and in the DDP
route (time ~-1.1, pH ~-0.1). Reading: a FAST 3-DX exit makes 3-DX quasi-steady and flattens the time course, but a
pH-FLAT exit then removes the ratio's pH dependence; both together need a fast 3-DX exit that grows with pH -- a
BASE-FAVOURED 3-DX sink (the hub review's H2, a cysteine trap of 3-DX with a free-amine factor, and signature vii's
pH 8 point). Expectations (sec. 18): joint pass in 600 draws (P 0.3) no; time-reaching draws with caramelisation
< -1 in >= 2 of 3 (P 0.6) no (median -0.53; one of three).

**v4 (arm B7ASST2M = B7ASST2 + r_tdp_mgo, the coded base-favoured 3-DX -> MGO exit; the engine's "base" factor, 1 at
pH 5).** Sampled additionally: log10 k_tdp_mgo (100 C, pH 5) in [-5, -1]; all other bands as v3; 600 draws, seed
20261102; same thresholds. This is a STAND-IN for any base-favoured 3-DX sink (a cysteine trap would have the same pH
direction); if v4 is feasible, the next step is to tell the candidates apart (MGO from 3-DX is measured in the Zhai
pots; a cysteine trap consumes cysteine).
**Expected (written now):** feasible P 0.4; the joint passes need log10 k_tdp_mgo > -3 (P 0.7 given feasible); MGO
rows overshoot at the joint passes (P 0.5 -- would favour a cysteine trap over the MGO route).

## 20. Feasibility v4 FEASIBLE; the MGO check rules out its stand-in; v5 = a cysteine trap of 3-DX (2026-10-31 ~11:20)

**v4 (B7ASST2M, 600 draws, seed 20261102; ..._feasibility_b7asst2m.md): FEASIBLE.** Time slope >= 5 / 8 by 73 draws
(v3: 3), pH segments >= 4 / 6 by 30, (vii) by 200; time + pH jointly by 7 -- several of them also pass (vii) and both
TDP declines (the first structure to pass every partition family). In most joint passes the split law sits near
"off" (pH_c 8.3 - 9.9): the ratio's pH dependence comes from 3-DX's REMOVAL, not from the enolisation split.
Expectations: feasible (P 0.4) yes; joint passes need log10 k_tdp_mgo > -3 (P 0.7) -- 3 of the first 5, the others
use a fast pH-flat exit with the split -- mixed.
**The pre-registered discriminating check (no fit):** at the joint passes the 100 calibrated Zhai MGO rows overshoot
by +0.49 to +1.28 dex (mean; rms 0.68 - 1.38) wherever the 3-DX -> MGO route carries the exit, against -0.04 at the
B7ASS reference; Cerny pot D 0.98 - 1.00 throughout. **So the base-favoured 3-DX sink must NOT make methylglyoxal**
(expected outcome, P 0.5, "would favour a cysteine trap"). Four independent lines now point at one mechanism -- a
CYSTEINE TRAP OF 3-DX, faster at high pH, giving a non-volatile adduct: the fast 3-DX exit (time course), its base law
(pH segments), signature (vii) (xylose competes for cysteine and rescues 3-DX at pH 8), and the cysteine review's
missing non-volatile cysteine sink (H3, P 0.7 as a needed second term). Wei 2022 shows cysteine trapping deoxyosones
(but reports 1-DX lowered before 3-DX -- mildly against, P ~0.3; wei2022_extraction.md).

**v5 (design, before code).** A new switch, `TDP_CYSTEINE_TRAP` (OFF: bit-identical): r_tdp_cys_trap, TDP + Cys ->
CBX + 7 FRAG_C + FRAG_N + FRAG_S (carbon, nitrogen, sulfur balanced; the cysteine amine declared destroyed, as
r_mgo_cys_strecker), second order, k(100 C, pH 5) x the engine's "base" factor, Ea 64.08 kJ/mol (declared, as the MGO
Strecker of cysteine). Arm B7ASST2C = B7ASST2 + the trap; feasibility sampling adds log10 k_trap in [-5, 0] L/(mmol
min). **Decided by:** feasibility on (i), (ii), (vii) as before, THEN at the joint passes the must-not-break checks
that v4 failed or did not read: the Zhai MGO rows (mean within +-0.3 dex of the reference), free cysteine and TTCA in
Zhai 2021 (the trap should LOWER the cysteine excess -- the cysteine review's test), Cao's blank pots (within 0.1 dex),
the thiol groups no worse than the reference.
**Expected (written now):** feasible P 0.5; at the joint passes MGO within +-0.3 dex P 0.7; free cysteine in Zhai 2021
moved toward the data by >= 0.3 dex P 0.4.

## 21. Held-out tests of the cysteine-trap hypothesis, declared 2026-10-31 ~11:45 BEFORE the v5 result is read

The partition signatures (i), (ii), (vi), (vii) have been used to design each feasibility structure (v1 -> v5); they
are now DESIGN TARGETS, not tests, for the trap. Credit for the trap (beyond "feasible") comes only from data not used
to design it:
1. **Another laboratory:** Hofmann 1999 (deoxyosone_hofmann1999 table; authentic quinoxaline standards) -- the ratio's
   plateau (~3) and its time shape; reported at the joint passes.
2. **The cysteine fate it predicts:** free cysteine and TTCA in Zhai 2021 (route table) must move TOWARD the data (a
   trap consumes cysteine); Cao 2017's blank pots must stay within 0.1 dex.
3. **Zhai 2023 JAFC's TTCA +- Cys double ratio** (CF-071; 0.53 - 2.6, the conflict with Zhai 2019) -- reported.
4. **What it must not make:** the calibrated Zhai MGO rows (within +-0.3 dex of the reference); the Davidek 2008
   acetic plateaus and pH 8 / 6 ratio; Cerny pot D no worse than the reference; the v2 thiol groups no worse.
5. **A power check of the feasibility method itself** before any credit: the tester must REJECT a deliberately wrong
   structure (R3En's box already never reaches (ii) -- recorded) and should be run on a synthetic case where the true
   structure is known (to be designed).
P(the trap is the real mechanism of the partition, written now) ~0.35 -- it was found by searching against these
signatures, and free cysteine (its co-reactant) is itself 2 - 40x too high in the engine, so the trap may be
compensating the cysteine error: the trap and the cysteine fate are to be tested TOGETHER.

**Amendment to sec. 21 after an independent red-team review (2026-10-31 ~12:00), STILL before the v5 result is read.**
(a) **Independence, honestly graded.** Item 2 (Zhai 2021 free Cys / TTCA) is NOT held out -- the missing cysteine sink
was one of the four lines that motivated the trap (sec. 20). Item 3 (CF-071) is one laboratory disagreeing with
itself and cannot discriminate. Item 1 (Hofmann 1999, xylose + alanine in phosphate, no cysteine) tests the 3-DX EXIT,
not the trap (the trap is inert there). (b) **Two held-out results already point AGAINST the trap:** Wei 2022 (another
laboratory) lowers 1-DX before 3-DX; and sec. 2's own argument -- a pure extra cysteine sink pushes TTCA lower, while
TTCA is already low (-0.14 to -1.41 dex). (c) **The partition signatures cannot tell candidates apart:** v4 showed a
different chemistry (3-DX -> MGO) passing every partition family; ANY fast 3-DX sink faster at high pH passes them, so
a "feasible" v5 is expected (P ~0.6) and low-information; the discrimination rests on the must-not-break rows.
(d) **Numeric pass thresholds for the must-not-break checks:** Zhai MGO mean within +-0.3 dex of the reference;
Davidek 2008 acetic plateaus and pH 8 / 6 ratio within their v2 tolerances; Cerny pot D no worse than the reference;
the v2 thiol groups no fewer than the reference's; Cao blanks within 0.1 dex; Zhai 2021 free Cys AND TTCA both
within 0.3 dex (must-pass, not "moves toward") together with signature (iv). (e) **Negative controls of equal
dimension, declared tests of the feasibility method (they must FAIL):** the same sink on DPO instead of TDP; the trap
pH-flat; the trap with its pH law inverted. **Positive control (must pass):** synthetic data from a known structure (a
v4 joint pass) with the fitted series noise. Feasibility is reported as a pass fraction with a confidence interval
across seeds, not mere existence. (f) **P(the trap is the real mechanism) revised to ~0.2**; the trap and the cysteine
fate are one structure for any further test; the partition is single-laboratory, provisional evidence (212 of ~276
deoxyosone rows are Zhai's, 3-DX / 1-DX against a diacetyl internal standard without authentic standards; P(a
substantial part of (ii) or (vii) is a method artefact) ~0.25) -- by roadmap row 3's rule its pH law is a declared band,
not a fitted coordinate, until another laboratory confirms it. **No trap fit follows v5.**

## 22. Feasibility v5 result (read 2026-10-31 ~12:15, after sec. 21 and its amendment were written)

**v5 (B7ASST2C, the cysteine trap, 600 draws, seed 20261103; ..._feasibility_b7asst2c.md): NOT feasible.** Time slope
>= 5 / 8 by 37 draws (best 6), pH segments >= 4 / 6 by 15 (best 6), (vii) by 101; time + pH jointly 0 (time + (vii)
15; pH + (vii) 6). Expectations: feasible P 0.5 (mine), ~0.6 (the red team) -- NO. Contrast v4 (the same sink written
as 3-DX -> MGO): feasible. The likeliest reason is the coupling sec. 21 flagged: the trap's rate is proportional to
free cysteine, which the engine has at the wrong level and with the wrong pH profile (excess largest at pH 5.5), so
the trap's EFFECTIVE pH law is distorted by the cysteine error. **The trap, coupled to today's cysteine chemistry, is
rejected without a fit;** it returns only as one structure with a corrected cysteine fate. The must-not-break checks
are moot (no joint pass). The switch stays in the engine, off.

## 23. A power test of the feasibility method (red-team item 3), designed before it runs (2026-10-31 ~12:25)

**Synthetic truth.** The first v4 joint pass (arm B7ASST2M; its constants and TDP exit as recorded in
..._feasibility_b7asst2m.json) generates synthetic TDP and DPO values on exactly the pots and times of the deoxyosone
table, with log-normal noise sigma = 0.15 dex per row (declared; the fitted series noise of these figures is 0.3 - 0.6
dex AFTER per-figure scales, so 0.15 on a scale-free ratio is a mild case). **Positive control:** the tester on arm
B7ASST2M against the synthetic data must find it feasible (expected P 0.8). **Negative control:** the tester on arm
B7ASST2 (the same structure WITHOUT the base-favoured 3-DX sink -- v3's arm) against the same synthetic data must find
it NOT feasible (expected P 0.7). 600 draws each, the bands of sec. 19, seeds 20261104 / 20261105. If the positive
control fails, "infeasible" verdicts carry little weight; if the negative control passes, "feasible" verdicts do.
