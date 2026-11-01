# Cluster 5 pre-registration: CF-076 / CF-011 re-run with Kallen-corrected bands and two added rivals (2026-11-01)

Written BEFORE any simulation in this folder. Prior knowledge declared:
- EX-054 (`docs/history/cys_exchange_2026-10-31/`): README, prereg, verdict, search results read (E4 fails only T3; T3
  and T2 fail the no-exchange control NC2 too; E1 core-feasible only at c_aa >= 20-30 /M).
- Today's dossiers read: kallen1971b, kallen1971 (its arithmetic 3 is wrong), hou2026, guan2012 (sec. 1), zhang2025b
  (sec. 3), zhai2019 / cao2017 / zhai2021 methods (PDF text: quench and storage sentences).
- **Zhang 2025b looked at before this file (arithmetic only, no model):** blank Q = TTCA/(Cys Xyl) rises 20 -> 72 ->
  ~18 600 L/mol from 5 to 50 min while TTCA is flat 17.7 -> 18.2 mM from 30 to 50 min and free Cys / Xyl fall to ~1 mM.
  Lysine arms (x1-x16): TTCA falls ~13 -> ~0.3 mM between 10 and 40 min (first order ~0.10-0.14 /min, nearly
  dose-independent). EX-054 never listed Zhang 2025b among its datasets (an omission under the "every dataset scored or
  excluded" rule); it enters here as test T8 (blank) and as a band source (lysine arms).

Folder: this file, `cys_exchange_c5.py` (copy of the EX-054 module, edited), `search_c5.py`, outputs `results_c5.json`,
`search_c5.json`, `cf011_probe.py` / `.json`, `verdict.md`.

## 1. Question and the rivals (neutral)

Why does one pH-flat mass-action exchange TTCA <-> Cys + Xyl (K = Cao's ~60 L/mol) fail to fit Zhai 2021's dilute
fed-TTCA pots together with the concentrated formation pots (Zhai 2019, Cao 2017), and now Zhang 2025b?

| id | structure | chemistry behind it |
|---|---|---|
| E0 | engine exchange, k_r(100 C) 0.0286 /min x[0.5, 2], Ea 64.1, K = K_ref, pH- and solute-flat | lumped fit |
| E1c | catalysed exchange, Kallen-corrected (sec. 2): general-acid catalysis of carbinolamine dehydration by ammonium (amino acids) and phosphate, SATURATING (ceiling 1 + rho), absent where attack limits (pH < ~5.5); label pH | Kallen 1971b Tables III, Fig. 7, Table I carboxylate null |
| E2c | in-situ pH: K(pH) by Kallen 1971 species bookkeeping (unchanged); k_r(pH) with Kallen's SHAPE (rises with pH to 5.6, flat to 6.6, falls above), replacing EX-054's acid limb a in [0, 1] | Kallen 1971b derived k_r(pH) |
| E3 | E0 + a second Cys sink | unchanged |
| E4c | E2c + E1c catalysis evaluated at the in-situ pH | |
| **E5** | **secondary-amine nucleophilic catalysis**: k_r = k0 (1 + c_N [TTCA + CAMA]), harmonically capped by the ring-opening rate k_open; K unchanged; label pH | Kallen 1971b Table IV (morpholine: catalyst iminium, then transimination by cysteine); TTCA ring N and the Cys-Amadori N are the only secondary amines in these pots |
| **E6** | **quench / storage operator**: plain pH-flat exchange (k0 band as E1c's) PLUS post-sampling formation Cys + Xyl -> TTCA at storage temperature, an integrated dose D (L/mmol) applied to every sample after the pot state; the pot's own initial TTCA comes from the rest of the t = 0 dose | Zhai 2019: "immediately cooled in an external ice bath and stored at low temperature (4 C) for HPLC detection"; 10 uL injected, no dilution stated; Cao: tap-water cooling, injected neat (standards to 0.5 M) |
| NC1 | instantaneous equilibrium (k_r 3 /min, K 0.06) | must fail T4 and T8b |
| NC2 | no exchange | must fail T5 / T7 / T8a |

## 2. Bands (all from data OTHER than the scored tests)

**Task 1, corrected exchange bands from Kallen 1971 + 1971b (cysteine + formaldehyde -> TC, 25 C, I 1.0 M):**
- TC hydrolysis k_r [Z, from measured k_f (1971b Table I / Fig. 6) / K_app (1971 Table IV)]: 1.4e-5 /min (pH 4.47),
  5.4e-5 (5.57), 5.6e-5 (6.56), 1.9e-5 (7.45); 1.8e-6 (2.84). Check: relaxation k_r(1 + 2[TC]/[Cys]eq) reproduces
  1971 Fig. 4 (t1/2 71 min at pH 5.39 vs read 1.2-1.4 h). k_f per TOTAL cysteine 0.6-4 /M/s (pH 5.6-7.5); the 1.1-2.9e3
  /M/s are per free-base amine.
- **k0(100 C) band: numerically unchanged [1e-5, 1] /min, basis rewritten.** TC's k_r x Arrhenius(25 -> 100 C,
  Ea 40-100: x25.6 - x3330) = 3.6e-4 - 0.047 /min (pH 4.5) and 1.4e-3 - 0.18 /min (pH 5.6). Transfer to TTCA needs a
  factor phi = k_r(TTCA)/k_r(TC) that no data on disk give: K(TTCA, 90 C) is ~1e-4 x K(TC, 25 C), and the split of that
  gap between k_f and k_r is unknown (phi 1 if all in k_f, 1e4 if all in k_r; < 1 possible for an alkyl-stabilised
  iminium). Declared phi in [0.03, 5] -> band edges 1e-5 (slow: pH 4.5, Ea 40, phi 0.03) and ~1 (fast: pH 5.6, Ea 100,
  phi 5). **Transfer assumptions:** (a) the release step of TTCA is the same C-N cleavage from the carbinolamine as
  TC's (N-H ring); (b) Arrhenius with one Ea over 25-140 C; (c) I = 1.0 M KCl vs low-ionic-strength pots ignored.
  Finding to state: at Ea 64 and phi = 1, TC itself would release ~1e-2 /min at 100 C, 10x above the dilute-pot ceiling,
  so the dilute pots need phi < ~0.1 or Ea near 40.
- **c_aa (ammonium general acid), pH-dependent:** c_aa(pH) = c6 x 10^(0.9 (pH - 6)) x g(pH); c6 in [0.01, 3] /M
  (log-uniform; Kallen 1971b Brønsted line slope 0.62 through H3O+ and water, amino-acid pKa 8.9-9.6: 0.1-0.3 /M at
  pH 6, 0.9-2.6 at 7.0, 2.3-6.4 at 7.4, +-1 log). g(pH) = 1 / (1 + 10^(2 (5.5 - pH))): attack limits below pH ~5.5 and is
  not buffer-catalysed (null for cyanoacetate, formate, acetate 0.05-0.4 M, pH 2.3-5.0).
- **c_P (phosphate):** initial slope 270 per M HPO4 2- at pH 6.1-6.6 (Fig. 7) -> band [90, 800] /M x g(pH).
- **Saturation:** F_cat = 1 + rho x / (rho + x), x = c_aa [AA] + c_P [HPO4 2-]; rho in [2.5, 4.0] (Fig. 7 plateaus
  x3.5 - x5 at 0.2-0.3 M, pH 6.12-6.62). EX-054's linear form had no ceiling.
- **E2c k_r(pH) shape:** log10 k_r = log10 k_r45 + s_lo (min(pH, 5.6) - 4.5) - s_hi max(pH - 6.6, 0); s_lo in [0, 0.85]
  (TC: +0.54 per unit from 2.8 to 5.6, +-0.3, floor 0 = flat), s_hi in [0, 0.8] (TC: -0.5 per unit 6.6 -> 7.45).
  K(pH) bookkeeping as EX-054 (pKa_T 5.6-6.6, dpKa_CT 1.6-2.6; Kallen 1971, unaffected by the error).
- **E5:** c_N in [0.5, 500] /M of (TTCA + CAMA), log-uniform. Basis: Kallen 1971b Table IV (10 C, formaldehyde 0.02-0.025
  M): catalysed / uncatalysed rate per total morpholine = 55-70 /M at low dose, pH 6.38 and 7.22 [Z: column x
  alpha_free-base Cys (Kallen 1971 microscopic pKa, +0.37 for 10 C) / footnote k0], levelling to ~19 /M at 0.3 M.
  Transfer x0.01-x8: TTCA ring N pKa ~6.2 vs morpholine 8.9 (more free base at pH 5-6, x~100; weaker nucleophile and a
  less stable iminium, x~0.001-0.1); CAMA pKa ~7.5-8.5. k_open(120 C) in [0.08, 0.16] /min: Zhang 2025b lysine arms
  (x1 - x16 lysine, dose-independent TTCA collapse = an amine capturing every opened ring; NOT the scored blank);
  Ea_open = the draw's Ea. Effective k_r = 1 / (1 / k_lin + 1 / k_open).
- **E6 storage dose:** D0 = 1.225e-3 L/mmol is the whole second-order dose that makes 14.5 mM TTCA from 82.6 + 165.2 mM
  (Zhai 2019 t = 0: x = a b (r - 1) / (b r - a), r = e^(D (b - a))). Jiangnan: storage share f in [0, 1] (D = f D0;
  the pot starts heating with the TTCA from (1 - f) D0). Xie lab (Cao, Zhang): D_x in [0, D0], declared (same
  chemistry, no t = 0 datum). Reverse reaction during storage ignored (25 C k_r ~5e-5 /min on the analogue). Dilute
  pots get the same Jiangnan dose (effect < 2 %). Assumption: every sample of a series stored alike.
- **Zhang 2025b inputs (T8):** Cys, Xyl forced from the rows (piecewise-linear; t = 0 is the 40 / 40 charge); CAMA
  forced, its rise drained from TTCA; total phosphate 0.2 M; label pH 6.5; in-situ pH (E2c / E4c) 6.5 - delta (1 -
  Cys/40), delta in [0, 0.6] (0.2 M phosphate capacity ~0.1 M/pH vs <= 40 mM H+ released); heating as a step at 120 C.
- Common bands, E0 / E3 bands, Zhai 2019 / 2021 pH laws, K_ref [0.045, 0.075] L/mmol: as EX-054.

## 3. Tests (EX-054 T1-T7 unchanged; T8 added)

| test | data | observable | sigma | lab |
|---|---|---|---|---|
| T1-T7 | as EX-054 prereg | | | Jiangnan T1-T6, Xie T7 |
| T8a | Zhang 2025b blank Fig. 2b | log10 TTCA(10 min) (formation) | 0.12 (come-up time not stated: a 2-3 min ramp in a 10 min integral; read 0.012) | Xie |
| T8b | same | log10 TTCA(50) / TTCA(30) (release gate: free Cys, Xyl ~1-4 mM) | 0.02 (two reads +-0.3 mM = 0.010; floor 0.01; forced-input error) | Xie |

**Power:** EX-054 method (N_score 600, N_truth 80, pairwise power >= 0.8 for at least one non-control pair, else the test
is dropped and reported).

## 4. Decision rule (unchanged from EX-054)

- **admitted, provisional:** joint feasibility > 0 (max |z| <= 2, sampled or by search) on the scored tests of BOTH labs
  and joint log-marginal within 3 nats of the best rival;
- **rejected:** search best max |z| > 3, failing test(s) powered against that rival;
- **inconclusive:** otherwise.
- **Primary set:** all scored tests T1-T8. Pre-registered sensitivities (declared now, informed by EX-054, so not post
  hoc for this run): S1 without T3 (T3 failed NC2 in EX-054); S2 = EX-054's T1-T7; S3 core T4-T8; S4 = K released:
  K_ref band widened to [0.045, 30] L/mmol (log-uniform) for E0-type slow exchange (E6 without dose and E1c), on T1-T8
  -- a hypothesis generator only, never an admission (K's band would then come from the scored tests).
- **A test that fails NC2 and every rival** is read as a shared-skeleton failure and reported, not used to rank.
- Controls must fail as designed (NC1: T4, T8b; NC2: T5 / T7 / T8a), else nothing is read.

## 5. Expectations (before any run)

| item | P |
|---|---|
| controls fail as designed | 0.95 |
| T8b powered (fast vs slow exchange) | 0.9 |
| T8 (a + b) fails every listed rival under K_ref [0.045, 0.075] (fast formation needs fast k_f; flat late TTCA needs slow k_r; together K >~ 2 L/mmol) | 0.65 |
| E1c rejected (c_aa now <= ~3 /M at pH 6-7, catalysis saturating, none at Cao's pH 4.2-4.7) | 0.85 |
| E2c rejected | 0.85 |
| E4c rejected | 0.6; inconclusive 0.35 |
| E5 admitted 0.1; inconclusive 0.3; rejected 0.6 | |
| E6 admitted 0.05; rejected 0.75 (dose <= D0 adds <= ~10 mM at Zhai 2019's early points and cannot make the 70 C plateau) | |
| no listed rival admitted on the primary set | 0.8 |
| S4 (K released) makes a slow exchange feasible on T4-T8 | 0.35 |

CF-011 probe (forced-input, sec. 6) expectations: Xie / Jiangnan apparent k_form ratio at matched temperature is >= x3
(0.75); the ratio is similar at 90 C (Cao) and 120 C (Zhang) within x2 (0.5); the storage operator moves Zhai 2019's
apparent k_form by < 30 % (0.9).

## 6. CF-011 probe 1 (forced-input version; the engine version stays undone)

Apparent first-order CAMA formation from TTCA, k_app = [Delta CAMA + k_cd integral CAMA dt] / integral TTCA dt over the
earliest window with rising CAMA, k_cd from the common band (de Roos, Ea 85.7) at both edges; per pot: Cao 90 C
(30-120 min), Zhai 2019 70 / 90 / 110 / 130 C (20-60 min; 130 C 0-20 min), Zhang 2025b 120 C (10-25 min). Zhai's own
Arrhenius line through 70-130 C is interpolated to 90 and 120 C. Phosphate-catalysis reading predicts a temperature-
independent Xie / Jiangnan ratio; a barrier reading predicts a ratio that changes with T. pH differs (Cao 4.2-4.7,
Zhang ~6-6.5, Zhai drifting) and is stated beside each ratio. The storage operator's bias on Zhai's TTCA (<= D0 x Cys x
Xyl) is propagated as an interval.

## 7. Addendum A (written 2026-11-01 while the main run was in progress, BEFORE reading results_c5.json or search_c5.json)

**A further unkeyed dataset found on disk: Zhai 2021 Fig. 7a** (13C5-TTCA formed from 10 mM 13C5-xylose + the Cys
released from 10 mM TTCA, 100 C, pH 5.5 / 6 / 7 / 8, 0-120 min; `data/articles/zhai2021.pdf` p. 10654). Isotope
labelling, evidence kind 1, on exactly the dilute re-formation channel; listed by neither EX-054 nor the dossier
tables. My reading (400 dpi crop `img/fig7a.png`, one reader, +-0.0005 mM):

| pH | 20 | 40 | 60 | 80 | 100 | 120 min |
|---|---|---|---|---|---|---|
| 5.5 | 0.0029 | 0.0032 | 0.0041 | 0.0052 | 0.0032 | 0.0029 |
| 6 | 0.0038 | 0.0042 | 0.0057 | 0.0059 | 0.0061 | 0.0055 |
| 7 | 0.0073 | 0.0137 | 0.0142 | 0.0078 | 0.000 | 0.000 |
| 8 | 0.0097 | 0.0198 | 0.0062 | 0.000 | 0.000 | 0.000 |

**Quantitation caveat (decides how it can be used):** "The concentration of 13C5-TTCA was obtained by subtracting the
content of unlabeled TTCA in the TTCA-[13C5-Xyl] reaction system from the content of TTCA in the TTCA-Xyl reaction
system" (Methods). A difference of two ~6-10 mM measurements in two pots has a resolution of ~0.1-0.4 mM (the Fig. 4b
replicate SDs), 10-100x the plotted values; the plotted +-0.0005 bars cannot be the difference's error. So:
- **S5a (face value, sensitivity only):** T9 = log10 13C5-TTCA at 20 and 40 min per pH (8 components), sigma 0.15 dex
  (read 0.02-0.07, label-pool dilution <= 0.05, floor). Predicted by adding labelled TTCA L to the +Xyl dilute pots:
  dL/dt = k_f Cys Xyl - (k_r + l_T) L, l_T = the forced TTCA's non-exchange loss per TTCA; all xylose treated as
  labelled (<= 10 % bias).
- **S5b (resolution-limited, one-sided):** predicted 13C5-TTCA <= 0.3 mM at every t <= 60 min (fail if above, z =
  (pred - 0.3) / 0.1).
- Both are applied to the sampled draws (joint with the primary set) and to each rival's best primary / S3 search
  point; neither enters the primary verdict. If the authors confirm a direct m/z 259 quantitation, S5a becomes primary.
- Expectation: at face value S5a rejects E0 / E3 / NC1 (P 0.9) and leaves the slow-dilute rivals (E5, E6, E4c)
  unaffected (P 0.7), i.e. it adds an isotope-grade line to the existing T4 gate; S5b is passed by everything except
  NC1 / E0 at the fast edge (P 0.7).
