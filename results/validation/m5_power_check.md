# M5 power check: sizing the designed deciders (synthetic data only)

> **Case:** [M5_buffer_fat_water_activity.md](../../docs/reference/cases/M5_buffer_fat_water_activity.md) sec. 5.3
> (coordinator rulings 2026-11-06). **Numbers:** [m5_power_check.json](m5_power_check.json), written by
> [scripts/diagnostics/m5_power_check.py](../../scripts/diagnostics/m5_power_check.py) (seed 20261106; one core,
> niced, threads pinned to 1; 4,332 s).
>
> **Rules kept.**
> - Synthetic data only. No measured value of any unit is read, no reserved or held-out data is touched, and the
>   script imports nothing from the engine or `data/`.
> - Every truth is drawn from the case's declared bands.
> - This is a power check, not a comparison: it decides nothing about any rival and spends no data.
> - The decision rule is the case's (sec. 4.0, ruling 1): a rival is **falsified only if p < 0.05 under BOTH forms**.
>   Form P is studentized per point. Form Q is pooled variance.
> - A pair is **separated** on a dataset simulated from rival A when rival B is falsified and A is not. One designed
>   dataset per simulation, so Holm has m = 1.
>
> **Limits.**
> - Each part is a reduced forced-input mini-model whose rival structures are exact. Real data carry misspecification
>   (M1's sink law, heat-up, matrix), so these powers are upper bounds on the design's discriminating power, not
>   forecasts.
> - Monte-Carlo counts per cell: B 24 (variants 16), F 16, W 200. A cell's MC SE is up to ~0.1 for B and F, and ~0.035
>   for W.
> - "Type I" below is the truth's own falsification rate. It includes fitter failures as well as the 5 % test size.
>   It is reported, not hidden.

## B. Buffer catalysis (the E5 row; CF-011 arm C; CF-094; CF-120)

**Mini-model.** The step rates are read directly (forced input) for four readouts:
- F, Amadori formation;
- D, Amadori decay;
- E, 2,3-enolisation;
- U, the furfural / FFT branch.

The design has five phosphate levels (inert zero, 5 mM, 50 mM, 0.2 M, 0.5 M) at two pH values (5.5, 7.0). Noise is
lognormal per replicate reaction, with CV 8 % (Zyzak's run-to-run ~6 %, range 2 - 12 %) or 15 %. Truths are drawn from
the bands:
- kappa per step U[5, 40] /M on [H2PO4-] at pKa2 7.0;
- R2c at K = 0.5 M with the same initial slope;
- R7: a jump U[0.5, 3] at K U[1, 5] mM, plus linear;
- R4': no factor on U.

Nested pairs are read one way only: R4' is inside R1, and R1 is inside R7.

**Full design, P(separated)** (CV 8 % / 15 %):

| pair (truth -> rival falsified) | n = 2 | n = 3 | n = 6 |
|---|---|---|---|
| R1 true -> R4' falsified (the step set) | 0.96 / 0.83 | 0.96 / 1.00 | 1.00 / 1.00 |
| R1 true -> R2c falsified (saturation at 0.5 M) | 0.88 / 0.29 | 0.96 / 0.88 | 1.00 / 1.00 |
| R2c true -> R1 falsified | 0.75 / 0.08 | 1.00 / 0.58 | 0.96 / 1.00 |
| R7 true -> R1 falsified (the 5 mM jump) | 1.00 / 0.96 | 1.00 / 1.00 | 1.00 / 1.00 |
| R1 ([H2PO4-]) true -> total-phosphate R1 falsified | 0.96 / 0.54 | 0.96 / 0.96 | 1.00 / 1.00 |
| total true -> species R1 falsified | 1.00 / 0.50 | 0.96 / 1.00 | 1.00 / 1.00 |
| type I (truth falsified), worst cell | 0.04 / 0.08 | 0.04 / 0.04 | 0.04 / 0.00 |

**Which arms are necessary** (n = 3, CV 8 %; each row drops one element):

| design | the pair that collapses | P(separated) without it |
|---|---|---|
| no 0.5 M arm | R1 vs R2c | 0.25 (R1 true) / 0.56 (R2c true) |
| no 5 mM arm and no inert zero | R7 vs R1 | 0.00 |
| one pH only (7.0) | species vs total | 0.00 |
| no furfural-branch readout | R1 vs R4' | 0.00 |

**Minimum size (B):**
- 10 conditions: 5 phosphate levels x 2 pH, plus the 0.2 M citrate arm, which serves CF-094's carboxylate question and
  is not sized here.
- All four readouts in the same pots.
- n = 3 replicate reactions per condition at a run-to-run CV <= 8 % (30 reactions): P >= 0.96 for every pair.
- n = 6 if the CV is ~15 % (60 reactions).
- Every arm is load-bearing: dropping any one collapses its pair.

## F. Lipid x thiol (the E7 row; CF-031)

**Mini-model.** Two phases (oil volume fraction 0.10). MFT and FFT pools are made at a forced rate and lost by:
- a water-phase sink;
- capture by enal (R2, R3);
- oil-phase oxidation (R6).

Partition is a fast equilibrium with K_MFT and K_FFT drawn from the declared log K band [1.0, 2.5].

The six arms are (a) buffer, fresh oil under air, fresh oil under N2; (b) enal dosed without oil; (c) PV-spiked oil
under N2, and enal dosed into oil under N2.

The readouts are MFT and FFT per phase, (d) thioether and oxidised products, and (e) a direct K measurement (CV 20 %).
The chemistry's fresh-oil fall is drawn as x3 - x15 (Yiltirak's x10 lies inside). Noise CV is 15 % or 25 % per cook.

R1 is nested in R3 and in R6, so "R1 true" rows read only R2's falsification.

**P(separated), n = 2 / 3 / 6, CV 15 %** (25 % in brackets where it differs by > 0.1):

| design | R2 vs R1 | R6 vs R1 | R3 vs R1 | R2 vs R6 | R2 vs R3 | R3 vs R6 | worst type I |
|---|---|---|---|---|---|---|---|
| full (a) - (e) | 1.00 / 0.88 / 0.94 | 1.00 / 1.00 / 1.00 | 0.94 / 1.00 / 0.94 | 1.00 / 0.88 / 0.94 | 0.94 / 0.88 / 0.94 | 0.94 / 1.00 / 0.94 | 0.12 (0.19) |
| no products (d) | 1.00 / 0.94 / 0.94 | 1.00 / 1.00 / 0.94 | 0.81 / 0.94 / 1.00 (0.75 / 0.88 / 0.88) | 1.00 / 0.94 / 0.94 | 0.94 / 0.94 / 0.94 | 0.81 / 1.00 / 0.94 | 0.06 (0.12) |
| no K (e), no products | 0.94 / 1.00 / 0.88 | 1.00 / 0.94 / 0.94 | 1.00 / 1.00 / 1.00 | 0.94 / 0.94 / 0.88 | 0.94 / 1.00 / 0.88 | 1.00 / 0.94 / 0.94 | 0.12 (0.19) |
| arms (a) + (b) only, no (c) - (e) | 1.00 / 0.81 / 0.81 | 0.94 / 1.00 / 0.94 | 0.81 / 0.94 / 0.88 | 1.00 / 0.81 / 0.81 | 0.88 / 0.81 / 0.81 | 0.81 / 0.94 / 0.88 (0.56 / 0.88 / 0.81) | 0.19 (0.25) |
| arm (a) only | 0.56 / 0.56 / 0.62 | 0.94 / 1.00 / 0.81 | 0.75 / 0.81 / 1.00 | 0.56 / 0.56 / 0.62 | 0.00 | 0.00 (R6 true); <= 0.25 (R3 true) | **0.44 - 0.69 (R2 true)** |

Each cell gives P(the first-named rival true -> the second falsified). Where both directions are meaningful
(non-nested), the lower of the two is shown.

**Reading (F).**
- **Arm (a) alone is not a usable design.** It cannot separate R3 from R6 (<= 0.25) or R2 from R3 (0.00), as the case
  states: the N2 arm separates the chemical routes from the physical ones, not from each other. In this mini-model it
  also leaves R2's dose and peroxide constants unidentified. The fitter then fails on R2's own data (type I
  0.44 - 0.69), which is a property of the ill-posed design, not of R2.
- **Arms (a) + (b) + (c), read per phase** separate every pair at P >= 0.88 at n = 3 cooks per arm (>= 0.81 at n = 2)
  at CV 15 %.
- The product readout (d) and the K measurement (e) add nothing in the exact model. They are kept in the design as
  protection against misspecification: they are the only readouts that test a route by its own product.
- **Minimum size (F):** six arms x n = 3 cooks (18 cooks), free + disulfide MFT and FFT per phase. At CV 25 % the
  powers stay >= 0.75 at n = 3.

## W. Water activity (CF-104's decider)

**Mini-model.** k_F of acrylamide formation at moistures measured directly (% wb, so no isotherm is needed for this
design). The sugar is 25 % of the dry matter. Noise is CV 20 % or 30 % per cook.

Rivals:
- R2: whole-sample (wet-basis) mass action, the shipped lane's basis;
- PW: per-water basis;
- R4: dissolved sugar capped at a solubility drawn U[0.15, 1.0] g / g water, so the cap binds inside the range;
- BAS: wet basis with a Bassama-like x2.6 slowdown at 70 % wb.

PW is nested in R4 and R2 in BAS, so those pairs are read one way only.

**P(separated), 5 moistures (10 - 70 % wb), CV 20 % (30 %):**

| pair | n = 2 | n = 3 | n = 4 | n = 6 |
|---|---|---|---|---|
| R2 vs PW (the basis; lower direction) | 0.81 (0.52) | 0.98 (0.94) | 0.98 (0.97) | 0.98 (0.97) |
| BAS true -> R2 falsified (is Bassama's wet-end slowdown real?) | 0.41 (0.15) | 0.89 (0.57) | 0.97 (0.92) | 0.96 (0.95) |
| R4 true -> R2 falsified (the solubility cap) | 0.26 (0.13) | 0.69 (0.38) | 0.80 (0.62) | 0.94 (0.82) |
| R2 true -> R4 falsified | 0.07 (0.03) | 0.26 (0.08) | 0.52 (0.20) | 0.82 (0.41) |
| type I, worst cell | 0.02 | 0.04 | 0.03 | 0.04 |

With 3 moistures (10, 40, 70 % wb) the basis and Bassama rows are within 0.05 of these. The R4 rows are 0.05 - 0.1
lower.

**Minimum size (W):**
- n = 3 cooks at each of >= 3 moistures (9 cooks) at CV <= 20 % decides the concentration basis (wet vs per water) and
  whether Bassama's wet-end slowdown is real (P >= 0.89). At CV 30 % this needs n = 4.
- Separating the solubility cap from the wet basis needs n = 6 at 5 moistures (30 cooks; P 0.82 - 0.94 at CV 20 %).
- The thiol a_w series has no written design and is not sized.

## Summary for the case (sec. 5.4)

| part | designed decider, minimum size | P(separated) at that size | what it cannot do |
|---|---|---|---|
| B | 10 conditions (5 phosphate levels incl. inert zero and 5 mM, x 2 pH) x 4 readouts, n = 3 reactions at CV <= 8 % (30), n = 6 at 15 % (60) | >= 0.96 every pair | transfer above 130 C; a ceiling in concentrated foods (needs C-B2's design) |
| F | 6 arms, phase-resolved MFT / FFT, n = 3 cooks (18) | >= 0.88 (>= 0.75 at CV 25 %) | misspecified sinks (M1) are outside the mini-model; arm (a) alone is not enough |
| W | >= 3 moistures x n = 3 cooks (9) for the basis and the wet-end slowdown; 5 x 6 (30) for the solubility cap | >= 0.89; cap 0.82 - 0.94 | the sulfur-lane a_w question (no design) |

**Re-check.** Every number here rests on the two-form rule of ruling 1. The M1.1 measurement-model revision (judge
calibration, in progress) must be applied, and this check re-run, before any of it freezes a design.
