# B65 plug-in mechanism-signature screen (lesson 2026-11-05)

*2026-11-25. Read-only diagnostic, run before `kinetic_core_b65_cysteine_redox_prereg.md` is rewritten. Engine as
shipped (`core_parameters("sulfur")`, B60 fit report), a handful of runs at shipped values and band corners, and
three per-pot or post-hoc plug-ins where the engine lacks the structure. Nothing fitted or sampled; no file under
`src/` or `data/` changed. Every pot is isothermal with the pH clamped at its stated value.*

- Script: `scripts/generators/screen_b65_cysteine_redox.py`
  (`./scripts/docker_maillard.sh run "python -m scripts.generators.screen_b65_cysteine_redox"`, 33 s on 8 workers).
- Numbers: `results/validation/kinetic_core_b65_plugin_screen.json` (all 88 variants x 57 pots, raw end states).

## 1. Pots and conventions

| pot | charge (mM) | T, time, pH |
|---|---|---|
| Zhang Fig. 1 | THI 44.5, PENT 99.9, OX 1.0 (ambient); + Cys 123.8 / GSH 48.8 / cystine as **Cys 124.8 + OX 62.4** (the b2_3 convention, `generate_kinetic_core_b2_3_fit.py`:412-416) | 115 °C, 60 min, 4.9 |
| Zhang Fig. 2d | Cys 4.1-619 mM (0.5-75 mg/mL); THI 10 / 30 / 100 / 133 mM (133 = Met 15 mg/mL at 1:3:3); xylose = thiamine by weight; **Met omitted** (no Met species in the sulfur lane) | 115 °C, 60 min, 4.9 |
| Grosch T. III | THI 0.060 (or 0.012) mM, ribose 1.34 mM, Cys 0-5.44 mM, OX 1.0 | 100 °C, 120 min, 5.7 |

Observables: free MFT; dimer share by mass MFTD/(MFT+MFTD). Measured: S1 ratios 19.1 / 15.0 / 28.7; S2 shares 0.36 /
0.079 / 0.063 / 0.35 (control / Cys / GSH / cystine); S3 MFT x27.5 from 0.5 to 62.5 mg/mL, flat to 75; S4 optimum
0.65 mM, max/Cys-free 3.7, 5.44 mM/Cys-free 0.38. Tolerances are the draft's (sec. 5). S4/S5 are spent: report only.

**Oxidant plumbing.** The engine has no constant supply term. "Supply s" is emulated with B11's own reservoir: OXR =
s x hold time, `k_ox_supply` = 10 (declared fast), so OX is held at OX(0) until the budget is spent. k_cys_ox is in
B11's units (L per ambient unit per min) with the engine's thiolate factor.

## 2. Results (selected rows; every row in the JSON)

S1 = MFT ratio to control, Cys / GSH / cystine. S2 = dimer share, control / Cys / GSH / cystine. S3 = MFT(62.5)/MFT(0.5)
at THI 10 / 30 / 100 / 133 mM. P = pass, F = fail, (p) = passes at some thiamine levels.

| # | variant | S1 | S1 | S2 | S2 | S3 | S3 | S4 (opt mM; 5.44/0) |
|---|---|---|---|---|---|---|---|---|
| 1 | shipped | 1.20 / 1.04 / 0.80 | F | .033 / .041 / .033 / .683 | F (control) | 1.26 / 1.30 / 1.24 / 1.21 | F | 5.44; 1.91 F |
| 2 | sugar route −1.9 dex | 1.00 / 1.00 / 0.70 | F | .033 / .033 / .033 / .637 | F | 1.00 (all) | F | 5.44; 1.01 F |
| 2 | sugar route +0.55 dex (Laplace 2σ of k_nf_mft) | 1.68 / 1.13 / 1.01 | F | .033 / .063 / .035 / .755 | F | 1.75 / 1.96 / 1.78 / 1.69 | F | 5.44; 4.2 F |
| 2 | sugar route +1.9 dex | 6.76 / 2.94 / 2.49 | F | .033 / .337 / .081 / .933 | F | 1.97 / 3.02 / 3.54 / 3.71 (step −0.15) | F | 5.44; 66 F |
| 2 | sugar route alone (k_hmp_mft = 0), shipped | (control 0) | – | – | – | 6.0 / 60 / 720 / 1200 | (p) | – |
| 3 | k_cys_ox 10^−5, closed or B11 reservoir | 1.20 / 1.04 / 0.80 | F | as shipped | F | as shipped | F | F |
| 3 | k_cys_ox 10^−3, closed | 1.21 / 1.04 / 1.02 | F | .033 / .003 / .033 / .296 | F (control) | 1.26-1.31 | F | F |
| 3 | k_cys_ox 10^−1, closed | 1.21 / 1.04 / 1.10 | F | .033 / .000 / .033 / .000 | F | 1.25-1.31 | F | F |
| 4 | OX(0) 0.01 or 0.1, any supply, k_cys_ox 0 | 1.20 / 1.04 / 0.79 | F | ≤ .004 / ≤ .004 / ≤ .003 / .68 | F | as shipped | F | F |
| 4 | OX(0) 10, any supply, k_cys_ox 0 | 1.18 / 1.03 / 0.83 | F | .247 / .293 / .251 / .707 | F (Cys, GSH) | 1.18-1.28 | F | F |
| 4 | OX(0) 10, supply 1, k_cys_ox 10^−1 | 1.14 / 1.03 / 0.85 | F | .247 / .000 / .251 / .522 | F (GSH only) | 1.31-1.45 | F | 5.44; 1.00 F |
| 4 | OX(0) 10, supply 0.5, k_cys_ox 10^−2.5 | 1.22 / 1.03 / 0.97 | F | .247 / .033 / .251 / .367 | F (GSH only) | 1.29-1.48 | F | F |
| 7 | OX(0) 10, no supply, k_cys_ox 10^−2.75 (GSH proxy .149) | 1.29 / 1.03 / 1.10 | F | .247 / .007 / .251 / .215 | P with the GSH proxy; F (engine GSH) | 1.30 / 1.38 / 1.48 / 1.50 | F | 5.44; 1.14 F |
| 4 | OX(0) 30, no supply, k_cys_ox 10^−2.5 | 1.41 / 1.03 / 1.18 | F | .479 / .006 / .484 / .266 | F (GSH only) | 1.37-1.88 | F | F |
| 5b | g(c), Kd 30 (any a, Kp corner) | 0.32-0.82 / 0.58-2.9 / 0.27-0.59 | F | Cys .013-.029 | F | 0.04-0.21 (falls) | F | 5.44; 2.5-10 F |
| 5b | g(c), Kd 0.3, a 1 | 0.21 / 0.04 / 0.18 | F | .033 / .010 / .001 / .365 | F (control) | 5.0 / 19 / 23 / 20 | (p)\* | 0.17; 0.91 P |
| 5b | g(c), a 10, Kd 1000 (what B1T needs at 115 °C) | 8.7 / 8.4-8.7 / 3.0 | F (Cys −0.34 dex) | .033 / .25 / .24 / .91 | F | 0.85-1.27 | F | 5.44; 7.7-10 F |
| 5c | slowed MFT loss, Kc 1 / 10 / 100 mM | Cys 5.2 / 4.5 / 2.2 | F | Cys .120 / .109 / .068 | F (control) | 1.7-2.0 / 3.2-4.4 / 3.1-3.9 | F | 5.44; 5.0 / 2.6 / 2.0 F |
| 6 | sugar +1.9 dex & Kc 10 | 20.9 / 7.8 / 3.0 | F (cystine) | .033 / .537 / .165 / .944 | F | 14.3 / 14.6 / 9.5 / 8.8 | (p) | 5.44; 92 F |
| 6 | sugar +1.9 & Kc 10 & OX 10, supply 1, k_cys_ox 10^−1 | 13.8 / 4.6 / 1.5 | F (GSH, cystine) | .247 / .000 / .602 / .708 | F (GSH) | 31 / 33 / 34 / 32 | **P** | 5.44; 1.6 F |

\* Passes S3 for the wrong reason: g at Kd 0.3 switches the thiamine route off (g ≈ 10^−2 to 10^−6 over the series) and
the rising MFT is the sugar route's; S1 is then −2 dex.

**S5** (12 µM thiamine, optimum near 0.68 mM): fails in every variant. Only g(c) at a 10, Kp 0.03, Kd 0.3 turns down,
and it peaks at the lowest dose.

**Shipped route attribution** (flux budget, Fig. 1): thiamine HMP → MFT 0.0286 mM in every arm. Sugar route adds
0.0070 (Cys, 20 %), 0.0008 (GSH, 3 %), 0.0070 (cystine), 0 (control). The engine's cysteine effect on MFT is
therefore the sugar route alone: x1.2.

### 2a. Variant 5a: fast exchange MFTD + 2 RSH ⇌ 2 MFT + RSSR at equilibrium, applied post hoc

K = [MFT]²[RSSR]/([MFTD][RSH]²) in {0.1, 1, 10}, on the shipped end states and on the OX 10 / supply 1 /
k_cys_ox 10^−1 end states.

- **Cys and GSH arms.** The share falls to 10^−8 to 10^−14 for every K. That passes S2's "< 0.15", but only because
  free thiol (77 mM) exceeds MFT (≈ 10^−5 mM) by 10^7. At equilibrium MFTD/MFT = [MFT][RSSR]/(K[RSH]²).
- **Cystine arm, engine convention** (78 mM free Cys left, RSSR = 62 mM): the share is 3 x 10^−3 to 3 x 10^−5, so it
  fails. It passes only if the dissolved free thiol is ≤ 3.2 / 1.0 / 0.24 mM (K = 0.1 / 1 / 10), i.e. cystine charged
  as disulfide with almost no free Cys. Then nothing in the arm makes MFT x29 (premortem cause 5).
- **Control.** There is no thiol partner, so exchange cannot set its 0.36. The oxidant must do that (OX(0) ~ 10 mM; row
  4).
- **The same equilibrium against the dimer shares on record** (my arithmetic, free thiol taken as charged):
  - Zhang Fig. 2f, Cys arm at 62.5 mg/mL (516 mM Cys, MFT 9.6 x 10^−5 mM, measured molar MFTD/MFT 0.031): needs RSSR ≈
    9 x 10^7·K mM. Impossible.
  - Xu 2010 (50 mM Cys, MFTD/MFT 0.40-0.58 by mass under N2) and Zhou 2023 (20 mM Cys, 6.5-9.6 %): with a generous MFT
    ≤ 10 µM and RSSR ≤ 10 mM, equilibrium allows ≤ 2 % at K = 0.1 and ≤ 0.2 % at K = 1.
  - **An exchange at equilibrium in the pot is refuted by every pot on record that has free cysteine and a measured
    dimer**, unless the dimer forms after the pot (measurement model M) or the exchange is slow (kinetic, not
    equilibrium-controlled).

## 3. Which structures move each signature the right way

- **S1 (MFT boost; x15-29).**
  - The thiamine route as shipped has no cysteine term. The x1.2 is sugar-route H2S.
  - Single levers, best corner:
    - sugar route +1.9 dex: x6.8 (+0.55 dex, its fitted 2σ edge: x1.7);
    - CF-089 g(c): x8.7, and only with Kd 1000 mM, i.e. no fall;
    - cysteine-slowed loss (C2): x5.2. Its ceiling is about k_mft_decay(115 °C) x 60 min = 0.143 x 60 ≈ 8.6, because
      a slowed sink only lets MFT approach its integrated formation.
    - k_cys_ox and the oxidant settings move S1 by < 0.1 dex.
  - **No single structure reaches ±0.3 dex on the Cys arm. Two do together** (sugar +1.9 dex x C2: Cys x20.9).
  - **The cystine arm fails in every variant** (best −0.97 dex). With cystine = Cys + OX, cystine can never beat
    cysteine, while measured cystine > cysteine. This is structural (charge convention and solubility), not a constant.
  - GSH reaches −0.17 dex only in the sugar +1.9 dex x C2 combination.
- **S2 (dimer share).**
  - Shipped passes the Cys, GSH and cystine arms. Cys and cystine are in-sample: the B2.3 / B60 fit rows are those
    ratios. It fails only the control (0.033 vs 0.36), and the Cys/control contrast goes the wrong way (x1.24 vs
    measured x0.22).
  - The control needs OX(0) ≈ 10 mM: 33x B11's saturation of 0.3 mM, with nothing on disk to support it.
  - **Cysteine scavenging (k_cys_ox ≥ 10^−3, i.e. inside B11's band) is what moves the Cys-arm share down. It passes
    control, Cys and cystine together** (rows 4, OX 10 / supply 0.5-1 or OX 30).
    - It fails GSH because the engine has no GSH + OX step.
    - With a GSH-as-Cys proxy (48.8 mM thiol), the four arms pass together at one point of a 27-point grid
      (OX(0) 5-20 mM, budget 0-0.25 mM/min, k_cys_ox 10^−2.75 to 10^−2.25; rows V7). That point is OX(0) 10, no supply,
      k_cys_ox 10^−2.75: control .247, Cys .007, cystine .215, proxy .149. It is a knife edge, with cystine and proxy
      each within 0.015 of their thresholds.
    - The reason is a thiol-equivalent budget. Cystine needs the oxidant to outlast its net 62 mM thiol, while GSH needs
      the oxidant to run out before 49 mM.
    - The GSH share is a figure read-off (`zhang2024b_extraction.md`: the text supports only the ordering GCys > Cys ≈
      GSH), so the edge is soft.
  - Equilibrium exchange (5a) moves Cys/GSH the right way but fails the cystine arm under the engine convention, and is
    refuted elsewhere (sec. 2a).
  - Sugar +1.9 dex and C2 move the Cys share **up** (dimerisation is second order in MFT): wrong direction.
- **S3 (x27 monotone rise).**
  - Shipped is flat (x1.2-1.3).
  - **The sugar route alone has the S3 shape and more** (thiamine route off: x6 to x1200, set by the xylose loading).
  - g(c) cannot rise more than x1.75 through the thiamine route (screen: ≤ x1.27 at Kd 1000; falls at Kd ≤ 30). The
    premortem's bound is confirmed.
  - C2 gives x2-4.4.
  - Only combinations pass: sugar +1.9 dex x C2 (x9-15), or that plus scavenging (x31-34 at all four thiamine levels).
  - **S3 therefore measures the sugar/thiamine split as much as any cysteine term on the thiamine route. It cannot
    discriminate Q-B unless the sugar route is pinned** (premortem cause 3: direction confirmed).
  - Model dimer-share direction over the series (measured Cys arm: 0.05 → 0.06, peak 0.15 at 50 mg/mL):
    - shipped: rises slightly (0.02 → 0.03, THI 30);
    - scavenging: falls to 0;
    - sugar +1.9 and C2: rise to 0.09-0.57, too high at the top.
- **S4/S5 (spent, report only).**
  - Every variant except g at Kd 0.3 rises monotonically with Cys in Grosch's pot. The sugar route +1.9 dex makes that
    x66.
  - g at Kd 0.3 gives the optimum and fall but costs −2 dex on S1.
  - Under every lever screened here, no constant set passes S4 with S1 or S3. The premortem's cause 2 holds: one
    T-independent concentration law cannot serve both regimes.
- **Measurement model (item 6, not modelled).** Zhang's read-out is 10 mL in a 40 mL vial, 45 °C, 30 min incubation
  plus 30 min SPME. Vessel and atmosphere are unstated for the reactor.
  - My arithmetic: 30 mL air at 45 °C holds 0.24 mmol O2, i.e. 24 mM O2 in the 10 mL sample (≈ 100 mM thiol-oxidising
    equivalents). That is far above the MFT (≈ 10^−5 mM), and of the order the in-pot model needs (OX ≈ 10 mM) to make
    the control's 0.36.
  - The record points the same way:
    - Xu 2010 finds MFTD/MFT 0.40-0.58 independent of N2/CO2/air and flags on-fibre oxidation;
    - Grosch's SIDA finds the dimer at ≤ 1-5 % of MFT;
    - Zhou 2023's share is pH-flat where in-pot thiolate oxidation would move ~100x.
  - **Vial-stage oxidation (with Cys scavenging in the vial) can plausibly set S2**, and nothing on disk separates it
    from in-pot chemistry.

## 4. Which structures survive the screen

| structure | verdict |
|---|---|
| A0-inert (shipped, k_cys_ox = 0) | Out on S1, S3 and the S2 control (needs OX ≈ 10 mM even then) |
| **A0 with k_cys_ox on its B11 band (scavenging)** | **Survives S2** for control, Cys and cystine with OX(0) ~ 10-30 mM. The engine's GSH arm fails for lack of a GSH + OX step. With GSH as a Cys proxy, all four arms pass only on a knife edge (OX 10, k_cys_ox 10^−2.75); soft, since GSH is a read-off. Moves nothing else |
| A1, exchange at equilibrium | Survives S2's Cys/GSH arms only. Cystine arm needs ≤ 0.24-3 mM free thiol. Refuted by Zhang Fig. 2f, Xu 2010 and Zhou 2023 dimer levels unless the dimer forms post-pot. Not worth sampling as an equilibrium; at most a kinetic exchange |
| B1 g(c), T-independent | Out: S1 ≤ −0.34 dex, S3 ≤ x1.27 via the thiamine route. Its only S3 pass switches the thiamine route off |
| B1T (Kd ≥ 1000 mM at 115 °C) | Out on S3 (≤ x1.27). Also needs E_d ≳ 280 kJ/mol (premortem) |
| C2, cysteine-slowed loss (premortem B5) | Moves S1 and S3 the right way, capped at about x8.6 by k_mft_decay·t. **Survives only combined** with a sugar route well above its fitted band |
| Sugar route on its fitted band (±0.55 dex) | Cannot move S1 or S3 enough |
| Sugar route at +1.9 dex | Needed for every S1 + S3 pass. Outside its Laplace band; its in-sample rows (Hofmann fed pots) were **not** re-scored here, so this corner may be excluded by the fit rows |
| Measurement model M | Not excluded by anything here. Plausible for S2 on the oxygen arithmetic |

**The combination that passes most:** sugar +1.9 dex x C2 (Kc 10) x scavenging (OX 10, supply 1, k_cys_ox 10^−1).
It passes S3 at every thiamine level, S1 Cys and S2 except GSH. It fails S1 GSH (−0.51), S1 cystine (−1.27) and S4.
It needs a sugar-route corner the fit rows probably forbid, and an oxidant 33x saturation.

## 5. What could not be run, and why

- No cystine (CSSC) or GSSG species, and no GSH + OX step in the engine. The cystine arm runs only under the b2_3
  convention; solubility-limited cystine appears only in the post-hoc 5a (RSSR = 2 mM).
- The exchange (5a), g(c) (5b) and C2 (5c) are plug-ins. 5a is post hoc on end states. 5b and 5c use the charged
  thiol as a constant per pot (cysteine loss over the run ignored; 38 % of the Fig. 1 Cys is gone by 60 min).
- Methionine in Fig. 2 is not modelled (no Met in the sulfur lane; it feeds MeSH/MMFT).
- B2 (HMP + OX), B3 (K_i on thiamine breakdown), A3 (abundance-weighted cross-coupling), the k_red_ox reductone buffer
  and the measurement model M were not run.
- No guard / in-sample fit rows were re-scored at any corner, so a "pass" here says nothing about breaking B60's rows.
- Zhou 2023 and Xu 2010 enter as arithmetic only (sec. 2a), not engine runs.

## 6. Pots seen by this screen (cannot be decisive later)

- Zhang 2024 Fig. 1a/b, all four arms. Already in-sample: the S2 ratios are B2.3 / B60 fit rows.
- **Zhang 2024 Fig. 2d (MFT) and 2f (MFT-MFT), Cys series, all ten doses**: run, and now SEEN. The GSH series of
  Fig. 2a-g was displayed while reading the dossier: SEEN, not run.
- Grosch 1992 Table III, both thiamine rows (already spent).
- Arithmetic only, also SEEN: Zhou 2023 dimer shares at pH 6/7/8 (already B17 T3), and Xu 2010 Table 1 MFTD/MFT.
