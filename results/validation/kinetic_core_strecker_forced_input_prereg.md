# Strecker aldehydes with FORCED INPUTS (roadmap step 2) -- pre-registration DRAFT (2026-10-31)

*Design only: no code, no fit, no engine call. Status: DRAFT, awaiting its independent review (rule: pre-registration
+ independent review before every fit). Written from `docs/guides/STATUS_AND_PLAN.md` secs. 4 - 5, the lane inventory
`kinetic_core_strecker_lane_prereg.md` secs. 1 - 7, the two-reader table `data/lit/fit_rows/strecker_hofmann2000.yml`
(111 rows) and its verification file, CF-020 / CF-070, the dossiers of Hofmann 2000 / 2000b, Cremer 2000, Balagiannis
2009, Chan 1994b, Weenen 2001, the ledger entries (trunk.yml, sulfur.yml), `FIT_HOLDOUT_DECLARATION.md` (no Strecker
amendment exists yet), the B29 pre-registration, and `src/kinetic_core/fitting.py`.*

**What this layer is.** Strecker aldehydes and acids predicted from MEASURED alpha-dicarbonyl courses, so the layer is
judged on its own: no deoxyosone hub (step 5), no composable network (step 6), no engine species. An error in the
dicarbonyl supply cannot propagate into it, and an error in it cannot be hidden by a supply constant. Scope: Phe ->
phenylacetaldehyde (PA) / phenylacetic acid (PAA) at 98 - 100 C, pH 7; Leu and Ala as identity tests. Out of scope:
pH law, temperature law (every row is 98 - 100 C), pyrazines from aminoketones (a sibling pre-registration).

## 1. The model

### 1.1 Rate law (one Strecker flux per donor, two exits)

For donor D in {MGO, GO, TDG (3-DG), G (glucosone)} and amino acid AA (Phe unless stated):

    J_D(t)      = k_D * f_AA * [D](t) * [AA](t)                          (the Strecker flux, L mmol^-1 min^-1)
    d[PA]/dt    = sum_D (1 - phi_D) * J_D  -  k_loss * [PA]
    d[PAA]/dt   = sum_D (phi_D + omega_D * O2) * J_D                      (O2 = 1 for "air", 0 for argon)

* **Acid not from the aldehyde.** PAA is a parallel exit, never an oxidation of PA: Hofmann 2000 Table 5 (labelled PA
  added to the glucose pot: < 1 % appears as labelled PAA in 30 min) and Table 4 (air raises PAA x1.6 - 5.5 while PA
  moves only x0.85 - 0.93). A PA -> PAA step is NOT in the model; a declared rival with it is not run (refuted by Table 5
  before any fit).
* **Oxygen adds acid flux, it does not divert the aldehyde.** If O2 split a shared intermediate, PA would fall by about
  the share PAA gains (MGO: x2 in total, PA should halve; it falls 8 %). So omega_D * O2 is an ADDITIONAL acid channel
  driven by the same donor x amino-acid encounter (e.g. oxidative decarboxylation of the transamination keto acid;
  mechanism not asserted). "air" is a categorical state of a closed vial with unstated headspace (O2 supply unknown, a
  dose-response is not available); argon is 0.
* **Amine form.** Total amino acid [AA]. Every fit and test pot is at label pH 7 (phosphate 0.1 - 0.5 mol/L), where the
  free-amine form f_NH2 = 1 / (1 + 10^(pKa - pH)) is a constant factor absorbed in k_D; the two forms differ only in
  pH extrapolation, which this layer does not do. Declared for later (step 3 / pH law): pKa(Phe alpha-NH3+) ~ 7.6 at
  98 C (9.13 at 25 C, ionisation enthalpy ~ 44 kJ/mol; approximate). Buffer molarity: Table 3 (0.5 mol/L, air) and
  Table 4 model B (0.1 mol/L, air) agree at 60 min within 0.02 dex (PA 1.06 / 1.02, PAA 2.06 / 2.04 mmol/L), so buffer
  catalysis of the Strecker step itself is below resolution; the SUPPLY's buffer catalysis is irrelevant here (forced).
* **Aldehyde loss.** First order, k_loss (aldol, imine / adduct formation lumped), one constant for PA; the acids are
  stable (declared). 3-methylbutanal takes PA's k_loss (declared +-0.5 dex); Balagiannis 2009's measured 3-MB sink
  (2.7e-3 /min at 130 C, Ea 78 +- 39) extrapolates to ~3.6e-4 /min at 98 C, a comparator for the fitted value, not an
  input. Over 30 - 60 min the loss is <= 3 %, so it matters only for Table 3's tail.

### 1.2 Fed pots (Tables 3, 4): the donor is NOT constant

The fed donor is charged at 100 mmol/L and never measured. It is consumed by the Strecker step and by every other
amine reaction (non-Strecker sink: imidazoles, melanoidins):

    d[D]/dt  = -kappa_D * [D] * [AA]        kappa_D >= k_D * f_AA;  s_D = k_D / kappa_D, the Strecker share
    d[AA]/dt = -nu * kappa_D * [D] * [AA]   nu = 1 declared (sensitivity nu = 0, AA constant)

(the A + B second-order solution is analytic; products by quadrature). Table 3's PA + PAA rises 1.63 -> 2.67 -> 3.12 ->
3.98 mmol/L (10 / 30 / 60 / 500 min): a plateau at ~4 % of the charge, so MGO's lifetime in this pot is ~15 - 25 min.
**Treating the fed pot as constant biases k_MGO low by x3 or more** (the 60-min chord 0.052 mmol/L/min vs the 10-min
chord 0.163); this, not only the yield-vs-rate confusion, is part of why B18's pyrazine-derived constant (log10 -7.53)
sits ~2.7 dex below the Strecker step. Rival donor-consumption law (reported, judged by `cost_drop_is_significant`):
amine-free first order, d[D]/dt = -kappa1_D [D]. The Table 3 Phe rows are NOT used to constrain kappa (sec. 3, CF-070).

### 1.3 Forced pots (the glucose + Phe pot)

[D](t) is the measured Fig. 3 course (TDG, MGO, GO, G at 10 - 500 min; [D](0) = 0), [AA](t) the measured Phe course
(hofmann2000b Table 3, amino-acid analyser, 100 C reflux: 92 / 82.5 / 64 / 51 / 39.8 mmol/L at 10 / 30 / 60 / 120 /
300 min), so PA and PAA are pure quadratures of the forcing -- no ODE, no supply constant, no thermal-history model (the
heat-up is inside the measured D(t)). Only k_D, phi_D, omega_D, f_AA and k_loss enter. Leu and Ala pots (Table 2) have no
measured donors: they are forced with the Phe pot's Fig. 3 courses (declared: amine identity changes the supply by
<= 0.3 dex, band carried into their prediction).

### 1.4 Two arms, both pre-registered, both reported

| arm | omega | fitted on | its oxygen hold-out | its glucose-pot test |
|---|---|---|---|---|
| **A** (the brief's split; minimal) | ONE omega for all donors | Table 3 (air, MGO) + Table 4 argon (4 donors) | Table 4 air GO / TDG / G | yes (secondary) |
| **B** (measured donor fingerprints) | omega_D per donor | Table 3 + Table 4 argon + Table 4 air | none (spent) | **yes: the decisive test** |

Arm A asks a mechanism question (is the oxygen acid channel donor-generic?). Arm B asks the layer's question (given
every fed-donor measurement, do measured dicarbonyls predict the sugar pot?). Neither is selected by its outcome: the
admission rule (sec. 8) is keyed to arm B now; arm A's verdict decides only the oxygen structure carried into step 7.

## 2. Coordinates and what the data can identify

| coordinate | arm | role | identified by | band |
|---|---|---|---|---|
| log10 k_MGO (98 C, pH 7, Phe) | A, B | FITTED | Table 3's rise (10 - 60 min) | [-7, -3] |
| log10 s_MGO (= k / kappa) | A, B | FITTED | Table 3's plateau | [-4, 0] |
| phi_MGO (argon acid fraction) | A, B | FITTED | Table 4 argon PAA / PA (+ Table 3 PAA) | [0, 0.9] |
| omega (A) / omega_MGO (B) | A, B | FITTED | Table 3 (air) vs Table 4 argon | [0, 5] |
| log10 k_loss (PA) | A, B | FITTED (weak) | Table 3's tail, 240 - 500 min | [-5, -2] /min |
| log10 s_D, D = GO, TDG, G | A, B | FITTED, exactly determined (one row each) | Table 4 argon PA + PAA at 60 min | [-4, 0] |
| phi_D, D = GO, TDG, G | A, B | FITTED, exactly determined | Table 4 argon PAA / PA | [0, 0.9] |
| omega_D, D = GO, TDG, G | B | FITTED, exactly determined | Table 4 air PAA | [0, 5] |
| **kappa_TDG** | A, B | **DECLARED** | Martins 2005 `k_tdg_mel` (3-DG + Gly, 8.12e-4 L/(mmol min) at 100 C, Ea 95.2, pH 6.8) x Gly -> Phe transfer | +-0.5 dex; hard floor below |
| **kappa_GO, kappa_G** | A, B | **DECLARED** | kappa_MGO (fitted) | +-0.5 dex; hard floor below |
| hard floors on kappa_D | -- | DERIVED from the forcing series | Fig. 3's post-peak net decline / mean [AA]: TDG >= 7e-5, GO >= 7e-5, G >= 1.7e-4 L/(mmol min) (consumption >= net decline; uses no aldehyde data) | -- |
| f_Phe | -- | 1 by definition | -- | -- |
| f_Leu, f_Ala | -- | **DECLARED 1** (identity-free null) | tested by Table 2 (sec. 3) | +-0.3 dex |
| nu (amine per donor consumed) | -- | DECLARED 1 | sensitivity 0 | -- |
| k_loss (3-MB) | -- | DECLARED = k_loss (PA) | Balagiannis comparator | +-0.5 dex |

**The minimum set the data identify.** From a single 60-min fed point, a donor's data give only its Strecker YIELD
(s_D, if the donor is used up) and its acid split -- not its RATE constant k_D = s_D * kappa_D, which is what a forced
pot at ~1 mmol/L donor needs. Only MGO has a fed time course (Table 3), so only k_MGO and kappa_MGO are separately
identified. For GO, TDG and G the rate rests on a declared kappa_D. **This is the layer's binding uncertainty**, and the
report must give its share of the glucose-pot predictive variance (sec. 6.4). Arm A: 11 fitted coordinates on 24 rows
(13 dof, all carried by the MGO rows: Table 3, Table 4 argon / air MGO); arm B: 14 on 30 (16 dof: those 13, plus one per other donor from its air PA). The exactly-determined coordinates earn no fit credit; they are
calibrations of the donors.

## 3. Rows: fit vs hold-out (declared before any integration)

**Blindness.** The inventory (sec. 3 - 4 of the lane pre-registration) and this design read Table 4's air / argon
ratios, Table 2's identity ratios and the pH shapes; sec. 7 below computes crude pre-fit numbers for the glucose pot.
Those hold-outs are therefore SEEN, not blind: the structure and the decision rule are fixed here, before any model
integration, but expectations are informed. The only rows nobody has read against a model are the ones still to be
keyed (hofmann2000 Fig. 2, hofmann2000b Fig. 1 / Tables 2 - 3, Cremer Fig. 6).

| rows | n | arm A | arm B | note |
|---|---|---|---|---|
| hofmann2000 Table 3 PA, PAA (MGO, air, 0 - 500 min) | 14 | FIT | FIT | SIDA |
| Table 3 PHE (semi-quantitative) | 7 | **excluded** | excluded | CF-070 duplicate candidate of hofmann2000b Table 3; diagnostic only (shape, tau = inf): the fitted nu * kappa_MGO predicts it |
| Table 4 model A (argon) PA, PAA, 4 donors | 8 | FIT | FIT | SIDA; buffer 0.1 mol/L as printed |
| Table 4 model B (air) MGO | 2 | FIT | FIT | near-replicate of Table 3 at 60 min (count once in any credit) |
| Table 4 model B (air) GO, TDG, G | 6 | **HOLD-OUT** (oxygen) | FIT | seen |
| Table 4 air + Cu (MGO) | 2 | reported | reported | no Cu term; direction only |
| Table 2 glucose + Phe PA, PAA, 30 min | 2 | HOLD-OUT | **HOLD-OUT: decisive (G1, G3)** | printed, SIDA; seen |
| Figs. 4 - 5 glucose + Phe pH 7 PA, PAA | 2 | HOLD-OUT | HOLD-OUT (G2, G4) | declared conditions (60 min, air): weight below G1 |
| Table 2 glucose + Leu (3-MB, 3-MB acid), + Ala (acetaldehyde) | 3 | HOLD-OUT (identity) | HOLD-OUT (identity) | forced by Fig. 3's Phe-pot donors (declared); seen |
| Table 2 acetic acid (Ala pot) | 1 | excluded | excluded | not Strecker-specific: the sugar makes acetic acid too |
| Fig. 3 dicarbonyls | 26 | FORCING (not scored) | FORCING | calibrated (quinoxalines) |
| Figs. 4 - 5 fed donors at pH 3 / 5 / 9 (+ pH 7) | 32 | signatures only (P group) | same | declared conditions; the layer has no pH law |
| hofmann2000 Fig. 2 (glucose pot PA, PAA 10 - 500 min; 6 x 2) | 12 | HOLD-OUT | HOLD-OUT (G5, G6) | **to key, two readers, before any fit** |
| hofmann2000b Table 3 (ARP-Phe, Phe in glucose pot) | 10 | forcing ([AA](t)) / diagnostic | same | to key; same laboratory |
| hofmann2000b Fig. 1 (glucose pot PA, PAA 0 - 500 min) | -- | HOLD-OUT, counted once with Fig. 2 | same | to key; "same pot in substance" (ledger) |
| hofmann2000b Table 2 glucose + Phe, argon / air / air + Cu, 120 min | 6 | PAA rows HOLD-OUT; PA rows **diagnostic only** | same | the PA air / argon ratios were FIT rows of B29 (refused): the never-promote rule bars them from credit |

**Cross-laboratory hold-outs.** None is keyed, and none on disk is a clean forced-input test:
* Pan 2025 (methional with MGO and 3-DG measured in the same pot, 100 - 140 C; figure-only, unread) is the only other
  pot with a measured donor and a Strecker aldehyde -- but its Table 2 constants were FIT in B22 (refused), so its pots
  can be scored only as a DIAGNOSTIC (100 C pots only, f_Met declared 1 +- 0.5 dex, methional's MeSH sink declared).
* Cremer 2000 (Ala / Val / Ile / Leu in ONE pot, aw 0.52, pH 5, headspace standard addition; Fig. 6 figure-only): the
  within-pot Ala : Leu aldehyde ratio is a pure identity ratio (one supply), a cross-laboratory signature of Hofmann
  Table 2's separate-pot Ala / Leu (0.54) -- jointly testing "supply is amine-independent". Needs keying (~2 h). Reported.
* Weenen 2001 (fed donors + Phe, pH 3.2, steam distillation; semi-quantitative): needs a pH law. **Candidate conflict
  (not entered; for the register):** Weenen's MGO / 3-DG phenylacetaldehyde ratio at pH 3.2 is 35, Hofmann's at pH 3 is
  ~1.6 (PA + PAA 1.47 vs 0.92 mmol/L); different MGO yields (21 % vs ~1.5 %), stripping vs closed vial.
* Balagiannis 2009, Huang 2017, Chan 1994b: no donor course (and Huang not simulable): comparators only.
* So the layer is, for now, a SINGLE-LABORATORY layer (DFA Garching). Gap G6 for the owner: pots that measure an
  alpha-dicarbonyl AND a Strecker aldehyde in the same pot, 80 - 130 C, water (no DOIs proposed: none verified).
  Gap G7: a fed GO / 3-DG / glucosone + amino acid pot with the DONOR followed over time (the kappa_D the layer
  declares); the cheapest designed experiment that would make the decisive test sharp.

## 4. The forcing and its error

* **Interpolation.** Primary: piecewise-linear in concentration through the measured points with [D](0) = 0.
  Alternatives scored as a discrete sensitivity: piecewise-linear in log10 concentration (from 10 min), monotone cubic
  (PCHIP). Glucosone's two zero-height markers (300, 500 min) are zero.
* **Error propagation (Monte Carlo, N = 500, seed declared).** Per draw: (i) each Fig. 3 point perturbed by its read
  error (0.001 - 0.005 dex; 0.29 / 0.44 dex for the near-axis glucosone points); (ii) one calibration scale per donor
  (multiplicative, whole series): N(0, 0.10 dex) for TDG / G (HPLC-DAD, synthesised reference quinoxalines, calibration
  not described), N(0, 0.15 dex) for GO / MGO (GC-MS against a labelled homologue, [13C4]butane-2,3-dione, no response
  factor printed); (iii) the run-to-run term between the Fig. 3 pot and the scored pot (Table 2 is a different run of
  the same recipe): the PA / PAA replicate floor (sec. 5); (iv) kappa_D uniform in its declared band (above its hard
  floor); (v) the fitted coordinates from their Laplace posterior. The 80 % and 95 % predictive bands come from the
  draws; the interpolation alternatives are reported beside them, not mixed in.
* **[AA](t).** hofmann2000b Table 3 (CF-070: if that series in truth belongs to the MGO pot, the glucose pot's Phe is
  unknown): sensitivity [AA] = 100 mmol/L constant. Effect on the 30-min prediction 0.06 dex, on the 120-min one <= 0.15
  dex (pre-fit arithmetic), so CF-070 is low-impact for the decisive rows.

## 5. Measurement model (`fitting.MeasurementModel`)

* PA, PAA, 3-MB, 3-MB acid, acetaldehyde: SIDA -> `scale_tau` 0 (absolute). One laboratory: `between_tau` empty (not
  estimable: one group); for any later cross-laboratory row the SIDA between-paper tau estimated in the other lanes is
  DECLARED, not re-estimated here.
* sigma per row = sqrt(read_error^2 + floor^2). The read errors (0.001 - 0.004 dex) are transcription only; no SD or n
  is printed anywhere. **Floor from the paper's own near-replicates** (Table 3 vs Table 4 air at 60 min; Figs. 4 - 5 pH 7
  vs Table 4 air, four donors; five pairs each): PA ~0.007 dex, PAA ~0.06 dex per measurement (pair RMS / sqrt 2). The PA
  pairs may be one experiment replotted (adjudication ruling 5), so the PA floor may be a re-reading floor: sensitivity
  with both floors x3.
* Model error on top: `fit_with_series_noise` with series POOLED by paper x species (`series_of = "hofmann2000|PHEAL"`
  etc.), so Table 4's single points share Table 3's noise estimate (a per-pot series of one or two rows gets none).
* The semi-quantitative Phe rows, if scored as a diagnostic: `semi_quantitative` (tau = inf, shape only).
* The two-level option is not needed for the fit (one absolute method class). It is the model for the reported Cremer
  identity ratios (outer group = paper x figure, `within_tau` per method), if Fig. 6 is keyed.
* Duplicates counted once (`count_once`): Table 2 Phe vs text p. 436; Table 5 PAA (0.295) vs Table 2 (0.29); Fig. 2 vs
  hofmann2000b Fig. 1; Table 4 air MGO vs Table 3 at 60 min in any credit count.

## 6. Signatures (frozen as v1 at review; hashed), chance baselines, nulls, power

### 6.1 The list (tolerance on |log10 predicted - log10 measured| unless stated)

| id | group | signature | tol | label |
|---|---|---|---|---|
| T1 | time | Table 3 PA + PAA shape, log10 P(t) / P(500), t = 10, 30, 60, 120 | 0.10 | in-sample |
| T2 | time | Table 3 PAA / PA, slope vs log10 t (1.63 at 10 min -> 2.26 at 500) | 0.10 per decade | in-sample |
| O1 | oxygen | PA(air) / PA(argon) for GO, TDG, G (measured 0.90, 0.93, 0.85) | 0.10 | held out (A); structural: passes by construction, chance ~1, NO credit |
| O2 | oxygen | PAA(air) / PAA(argon) for GO, TDG, G (measured 5.5, 2.2, 1.6) | 0.20 | held out (A) |
| O3 | oxygen | rank order of the O2 effect over the four donors (measured GO > MGO > TDG > G) | exact | held out (A) |
| G1 | glucose pot | PA + PAA at 30 min, Table 2 (0.071 mmol/L) | decision rule, sec. 8 | **held out, decisive** |
| G2 | glucose pot | PA + PAA at 60 min, Figs. 4 - 5 pH 7 (0.141) | x2 | held out (declared conditions) |
| G3 | glucose pot | PA / PAA at 30 min (1.45) | 0.15 | held out, scale-free |
| G4 | glucose pot | PA / PAA at 60 min (1.23) | 0.15 | held out, scale-free |
| G5 | glucose pot | shape log10 P(t) / P(120), t = 10 ... 500, Fig. 2 (to key) | 0.15 | held out, scale-free |
| G6 | glucose pot | PAA / PA at 500 min (text 1.7, hofmann2000b) and the PA = PAA crossover time (text: after ~60 min) | 0.15; x2 in time | held out, scale-free |
| G7 | glucose pot | the gap R(t) = measured / predicted, its trend 30 -> 60 -> 120 min (falls if the missing flux tracks ARP-Phe, which peaks at 30 min) | reported | diagnostic |
| I1 | identity | Leu / Phe (aldehyde + acid) at 30 min, Table 2 (1.94) under f = 1 | 0.30 | held out; a test of a declared value (no chance rate) |
| I2 | identity | Ala / Phe (aldehyde only) at 30 min (1.14) under f = 1 | 0.30 | as I1; acetaldehyde volatility (bp 20 C) noted |
| P1 - P4 | pH | PA peaks at pH 5 (all 5 systems); PAA / PA falls monotonically with pH; glucose pH 5 / 9 PA ratio (2.7) vs fed (1.8 - 2.3); fed totals flat pH 5 - 9 and lower at pH 3 (x1.2 - 2.1) | -- | recorded, NOT scored: the layer has no pH law; targets for step 3 |

### 6.2 Chance baseline
`--chance 200` in the style of `score_deoxyosone_signatures.py`: every coordinate uniform in its band (fitted ones in
their fit band, declared ones in their declared band), the forcing at its central reading; pass rate per signature and
per group. A group earns credit only if its chance rate is < 0.2 (O1 earns none by construction).

### 6.3 Nulls (credit requires beating both)
(a) The no-change / constant-ratio null of `build_contrast_scores.py`: a scale-free ratio predicted equal to the fed
pots' mean (PA / PAA = the mean of Table 4's donors under the same atmosphere; O2 ratio = MGO's). (b) The typical-yield
null for G1 / G2: the glucose pot's yield predicted as the fed pots' mean Strecker yield per mmol donor times the
time-integrated measured donor -- i.e. forcing WITHOUT rate constants (k_D replaced by the cross-donor mean yield /
(D0 * 60 min)); the layer's kinetics must improve on it to be worth having.

### 6.4 Power and variance decomposition
With one laboratory and ~6 credited held-out signatures in 3 groups, power is low; it is REPORTED as the share of
synthetic data sets (the fitted model + the estimated noise, 200 draws) on which each held-out group passes (the
"does the scorer pass the true model" rate), and the false-pass rate from the chance draws. A study bootstrap is
degenerate (one study): pots are bootstrapped within the study and labelled so. For G1 the report gives the predictive
variance shares of kappa_TDG, kappa_GO, kappa_G, k_MGO, donor calibration, interpolation and replicate floor.

## 7. Pre-fit arithmetic (scratch, no fit; disclosed because it shaped sec. 8's probabilities)

Trapezoid integrals of Fig. 3 times [AA](t), 0 - 30 min (mmol^2 L^-2 min): TDG 1591, GO 692, MGO 168, G 45; 0 - 60 min:
4673, 2278, 686, 151. k_MGO (air total) ~ 2e-5 L/(mmol min) from Table 3's first chord (a lower bound; B18: 3e-8).
* G1 under arm B: k_D from the Table 4 air yields with kappa_D = kappa_MGO gives PA + PAA ~0.024 mmol/L at 30 min
  (measured 0.071, R ~3) and ~0.08 at 60 min (measured 0.141, R ~1.8); the "no donor decay" limit (kappa_D -> 0) gives
  0.006 / 0.020; Martins' kappa_TDG (~1.4 kappa_MGO) moves the 30-min value to ~0.029; kappa_GO at +0.5 dex to ~0.048.
  So the central expectation UNDER-predicts by x2 - 3, and the declared kappa band alone spans about the measured value.
* G3 / G4 / G6 under arm B: PA / PAA ~1.17 (30 min; measured 1.45), ~1.10 (60 min; 1.23), ~0.7 - 0.8 (500 min; 0.59):
  the donor mix shifting toward MGO reproduces the acid's overtaking in direction.
* Arm A: omega ~1.0 (per unit argon Strecker flux, from MGO) predicts PAA(air) / PAA(argon) of x6 (GO, measured 5.5:
  pass), x7.6 (TDG, measured 2.2) and x8 (G, measured 1.6): O2 fails for both C6 osones and O3's order is inverted
  (the model ranks by 1 / phi_D). The same omega gives the glucose pot PA / PAA ~0.7 (measured 1.45; G3 fails by 0.3 dex).
  The oxygen acid channel is, on these numbers, specific to the C2 / C3 fragments -- arm A is a near-known failure whose
  value is to fix that structure for step 7.

## 8. Decision rule and expectations (probabilities written before any fit)

**Pre-conditions (no verdict without them):** two starts agree (cost within 1 %, identified constants within their
Laplace sigma); no budget-exhausted fit; T1 / T2 in-sample within tolerance after the noise estimate (else the fed-pot
structure is wrong and nothing downstream is read).

**The layer's verdict (arm B, on G1 with its 80 % predictive band; R = measured / predicted):**

| verdict | condition | P |
|---|---|---|
| ADMIT | R in [1/2, 2], and G3 passes | 0.30 |
| ADMIT AS LOWER BOUND | R in (2, 10] and the band's upper edge below the measured value: measured dicarbonyls do not suffice; the missing flux (ARP-direct oxidative route, hofmann2000b; unmeasured donors 1-DG / C4 - C5 fragments / DHA; donor speciation at 100 mmol/L) is quantified and becomes step 7's next arm | 0.25 |
| INCONCLUSIVE | R in (2, 10] but the band covers the measured value: the declared kappa_D band binds; remedy G7 (or a literature kappa), no re-fit to rescue it | 0.25 |
| REFUSE | R < 1/2 beyond the band (over-prediction cannot be fixed by adding routes: the rate does not transfer from 100 mmol/L fed donors to ~1 mmol/L, e.g. a donor-second-order path), or R > 10 | 0.20 |

Other expectations: P(converged) 0.9; P(T1 / T2 pass in-sample) 0.6; P(second-order donor consumption beats the
first-order rival by the chi-square rule) 0.5; P(log10 k_MGO in [-5.0, -4.4]) 0.7; P(k_MGO >= 2 dex above B18) 0.9;
arm A: P(O2 passes for GO) 0.7, for TDG 0.05, for G 0.05, P(O3 exact) 0.02; arm B: P(G3) 0.55, P(G4) 0.6, P(G5) 0.4,
P(G6) 0.45; P(I1) 0.55 (measured 0.29 dex against 0.30), P(I2) 0.75; P(G2 within x2) 0.5; P(R30 > R60, i.e. the gap
shrinks as ARP-Phe falls) 0.6. Any credited group must also beat both nulls (sec. 6.3).

## 9. What it needs in code (hours of work; compute is minutes)

A standalone module, not the engine: the engine has no Phe / PA / PAA species and the network refactor (step 6) is
pending; the forced pot is a quadrature and the fed pot an analytic A + B solution, so no stiff integrator is needed.
1. `fit_rows.load(..., extra_species=...)`: the Strecker layer declares its provisional keys (PHE, PHEAL, PHEAC, LEU,
   MEBUAL3, MEBUAC3, ALA, ACAL) instead of the loader refusing them; a unit test that the engine vocabulary is
   unchanged. ~1 h.
2. `src/kinetic_core/strecker_forced.py`: rate law, fed-pot solution, forced quadrature, interpolants, the Monte Carlo
   forcing; tests: quadrature vs analytic on a fed pot, "forced = fed" when D(t) is the fed solution, oxygen off = argon.
   ~4 h.
3. `scripts/generators/fit_strecker_forced_input.py`: arms A / B via `fitting.Problem`, `two_starts`, `laplace`,
   `fit_with_series_noise` (series pooled by paper x species), `count_once`; JSON with the arm read from the record
   (rule e). ~2 h.
4. `scripts/generators/score_strecker_signatures.py` (reusing `score_mechanism_signatures` helpers): v1 list, chance
   baseline, nulls, synthetic power, variance decomposition. ~3 h.
5. Keying by two readers + adjudication BEFORE the fit: hofmann2000 Fig. 2; hofmann2000b Tables 2 - 3 and Fig. 1
   (~3 - 4 h of agent time); Cremer Fig. 6 optional (~2 h). A FIT_HOLDOUT_DECLARATION amendment for sec. 3 (~0.5 h).
6. Independent review of this draft (~1 h); report, experiment-register entry, lessons (~1.5 h).
**Total ~16 - 19 h** (the roadmap's 17 - 21 h includes the pyrazine sibling); P(done by 11-07) ~0.65.

## 10. How this could mislead, and the remedies

1. **The forcing has its own calibration.** Fig. 3's dicarbonyls are quinoxalines (HPLC-DAD with synthesised references;
   GC-MS against a labelled homologue, no response factor printed), the responses SIDA: a donor-scale bias maps 1 : 1
   onto G1. Remedy: per-donor calibration draws (sec. 4); the scale-free G3 - G6 carry the decision's second leg; the
   variance decomposition says how much of R is calibration.
2. **The derivatisation measures "free + reversibly bound" donor after cooling** (OPD 3 h at 30 C, no trap in the pot):
   donor held as Phe hemiaminals at 98 C may be released (over-reading the reactive pool) or not (under-reading).
   Remedy: none in the data; noted as part of the calibration term; the ARP's breakdown during derivatisation is
   negligible (3 h at 30 C ~ 0.2 min at 98 C for a barrier near 100 kJ/mol).
3. **Unmeasured donors.** 1-DG, 3-deoxyerythrosone / -xylosone, DHA, glyceraldehyde, diacetyl are Strecker donors
   (Weenen's ranking) and are not in Fig. 3: forcing four donors is a LOWER bound on the dicarbonyl route. This is why
   under-prediction is not refutation and over-prediction is (sec. 8).
4. **A dicarbonyl-free route.** hofmann2000b: the Amadori compound gives PA oxidatively without a free dicarbonyl
   (air / argon x9 from ARP-Phe, x3.5 from glucose + Phe). In air the glucose pot carries it; the forced layer does not.
   Remedy: G7's trend against ARP-Phe's course; the argon glucose pot (hofmann2000b) as a diagnostic upper-bound check
   (its donors were not measured under argon, and its PA ratios were B29 fit rows).
5. **Yield is not rate.** Fed 60-min points measure Strecker SHARES; the forced pot needs RATES. The declared kappa_D
   bridges them and can dominate the prediction. Remedy: report the share; INCONCLUSIVE is a pre-declared outcome, not a
   licence to re-fit; G7 (fed donor followed) is the fix.
6. **Fed-donor speciation.** At 100 mmol/L, MGO / GO are hydrated and oligomeric; at ~1 mmol/L mostly monomeric
   hydrates: an apparent k per total donor need not transfer. Over- or under-prediction both possible. Remedy: none on
   disk (no donor-concentration ladder); flagged in any ADMIT.
7. **Seen hold-outs.** Table 4 air, Table 2 and the pH shapes were read before this design (sec. 3); sec. 7 computed
   rough answers. Remedy: credit only from groups whose chance rate is < 0.2 and that beat both nulls; the to-be-keyed
   rows (Fig. 2, hofmann2000b Fig. 1) are keyed blind to any model output and carry G5 / G6.
8. **Declared conditions.** Figs. 4 - 5 (60 min, air) are the adjudicator's inference; G2 / G4 are weighted below G1 and
   the P group is not scored.
9. **One laboratory.** Every fit and test row is DFA Garching: a shared analytical chain cannot be separated from the
   chemistry; leave-one-lab-out is impossible. Remedy: say so on every output; Pan 2025 / Cremer as diagnostics;
   G6 / G7 to the owner; no product-level number from this layer.
10. **CF-070 and the Phe course.** Low impact on G1 / G2 (<= 0.06 dex); the Table 3 Phe rows stay out.
11. **Ala's acetic acid and acetaldehyde volatility.** Acetic acid excluded (sugar origin); acetaldehyde (bp 20 C) in a
    refluxing pot may be lost to the condenser: I2's tolerance stands, a miss is read as volatility first.
12. **Arm forking.** Two arms could be read as choose-the-winner. Remedy: both are scored on every signature, the
    admission is keyed to arm B now, arm A decides only the oxygen structure.

## 11. For the reviewer and the owner

(a) Is arm B (Table 4 air spent) the right decisive configuration, given the brief proposed holding Table 4 air out?
(b) Is the declared kappa_TDG from Martins (glycine, pH 6.8) acceptable, or should all three kappa_D share one band?
(c) Enter the Weenen-vs-Hofmann pH-3 donor-ratio contradiction in the conflicts register? (d) Key hofmann2000 Fig. 2
and hofmann2000b before the fit (recommended: they carry the only blind scale-free tests). (e) The fit_rows loader
change (sec. 9.1) touches `src/`; it is the only engine-side edit.
