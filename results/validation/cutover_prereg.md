# Cutover final exam — PRE-REGISTRATION

**Written 2026-08-29, on `audit-remediation` @ `5a9135d`, BEFORE the scorer was written and
BEFORE any measured value in `data/benchmarks/external_validation/` was read.**

This document states, per bundle family, what the kinetic core is expected to do on the 21
frozen external-validation bundles. It is written from the four module scorecards
(`results/validation/kinetic_core_b{1,2_1,3,4}_holdout_report.md`) and from the core's own
species and reaction inventory — **not** from the bundles' measured columns.

What WAS read from the bundle files to write this: `benchmark_id`, `conditions`,
`precursors`, and the KEYS of `holdout_targets` (compound names only). No value under
`holdout_targets` or `reference_volatiles` was read, printed, or used. That is the minimum
needed to decide which lane a bundle routes to and whether its target compounds are core
species at all, and deciding that is the whole point of an envelope declaration.

Wave B5 is a WIRING AND SCORING wave. **No parameter changes in this wave, at all.** The
frozen values are exactly those in the B1/B2.1/B3 fit reports.

---

## 1. What the core actually is, and why that decides the envelope

The kinetic core is **three networks, not one**, and they do not compose:

| lane | network | integrator | species added | pH | scorecard |
|---|---|---|---|---|---|
| **L1 trunk** | `REACTIONS` (15 steps) | `integrate` | Glc, Fru, Gly, SB, AMA, TDG, ODG, MGO, FA, AA, MEL | **none** (pinned, `NETWORK_PH`) | B1: browning median **1.45x**, PASS |
| **L2 sulfur** | `FULL_REACTIONS` = trunk + sulfur (79 steps) | `integrate_sulfur` | PENT, Cys, THI, H2S, ARP, DPO, TDP, DDP, NF, FUR, MFT, FFT, MFTD, ACTZ, … | **yes**, incl. a pH trajectory | B2.1: **15/27** gating rows |
| **L3 acrylamide** | `FULL_ACRYLAMIDE_REACTIONS` = trunk + acrylamide (31 steps) | `integrate_acrylamide` | Asn, SBA, INT1, Asp, ACR, ACRCYS, Gln, Lys, Ala | **none**; a_w is metadata only | B3: **0/4** gating rows |

The two facts that do most of the work below:

1. **L2 and L3 are disjoint in their STEPS.** `FULL_ACRYLAMIDE_REACTIONS` deliberately omits
   every sulfur step (`acrylamide.OUT_OF_SCOPE`: composing them would spend the same cysteine
   twice). So no single run can carry both a sulfur product and acrylamide. A bundle whose
   targets span both lanes cannot be answered by one integration, and the engine must say so
   rather than pick a lane silently.
2. **The core has NO lipid-oxidation path and no hydroxymethylfurfural.** Verified by
   inventory, not assumed: `grep -c HMF src/kinetic_core/*.py` returns zero everywhere, and
   there is no hexanal, nonanal, 1-hexanol, 2-pentylfuran or DMHF species in any of the three
   state vectors. These are not "hard cases"; they are compounds the model cannot name.

## 2. Envelope declarations, decided in advance

A point is **OUT-OF-ENVELOPE** when its target compound is not a species in any core lane, or
when its precursors cannot be represented in the lane its target requires. An
out-of-envelope point is **declared, not predicted**. It is scored as neither a pass nor a
fail, and the report must not emit a number for it.

| bundle family | n bundles | n points | target compounds | lane | **declaration** |
|---|---:|---:|---|---|---|
| Hofmann 1998 cys/sugar 145 °C | 5 | 10 | FFT, MFT | L2 | **IN** |
| Yiltirak 2026 rib/cys 100–130 °C | 4 | 8 | MFT, FFT | L2 | **IN** |
| Chang 2021 / Ye 2024 glc+asn 180 °C | 3 | 3 | Acrylamide | L3 | **IN** |
| Lin 2022 fru+asn 180 °C | 1 | 2 | Acrylamide / **HMF** | L3 | **1 IN, 1 OUT** |
| Chang 2021 30 min water | (above) | 2 | Acrylamide / **HMF** | L3 | **1 IN, 1 OUT** |
| Schibilsky 2019 glc+ala 130 °C pH 5 / pH 8 | 2 | 6 | Furfural, DMHF, HMF | — | **OUT (all 6)** |
| Steinhagen 2021 glucose-only autoclave | 1 | 1 | HMF | — | **OUT** |
| Matrix path (Bi ×2, Liu, Li) | 4 | 8 | hexanal, nonanal, 1-hexanol, 2-pentylfuran | — | **OUT (all 8)** |
| **total** | **21** | **40** | | | **23 IN / 17 OUT** |

Named reasons for every OUT:

- **HMF (5 points: Lin, Chang-water, Schibilsky ×2, Steinhagen).** 5-hydroxymethylfurfural is
  not a species in any of the three lanes. The hexose dehydration route that makes it was
  never parameterised, because no dataset in the fit corpus measures it. Declared out; no
  number emitted.
- **DMHF (2 points, Schibilsky ×2).** 2,5-dimethyl-4-hydroxy-3(2H)-furanone is not a core
  species. The core carries **norfuraneol** (`NF`, 4-hydroxy-5-**methyl**-3(2H)-furanone),
  which is a *different compound* on a different (pentose) route. Predicting DMHF from the NF
  node would be a species substitution the corpus does not license. Declared out.
- **Furfural in the Schibilsky bundles (2 points).** `FUR` *is* a core species and
  `r_glc_fur` *does* make it from glucose without an amine — so a number is computable. It is
  nonetheless declared **out**, on the AMINE axis: Schibilsky's system is glucose +
  **alanine**, and the L2 network has no alanine species at all (`Ala` exists only in L3, as
  a competitor that consumes glucose). Substituting the trunk's glycine for alanine is not a
  free move — `species.py` warns in terms that an α-amino acid of the same C/N count
  "substitutes without changing the bookkeeping, but NOT without changing the rates", and
  alanine (C3N1) does not even match glycine's (C2N1) bookkeeping. The amine-free
  caramelisation number is reported as a **declared diagnostic**, clearly marked, and is not
  scored.
- **The 8 matrix-path points.** The core has no lipid-oxidation path. This is the case the
  build brief names explicitly, and the pre-registered answer is to declare them
  out-of-envelope rather than predict garbage. **This is also the honest reading of the old
  lane's 0/8**: the old lane produced eight numbers here; the core produces none, and
  produces none for a stated structural reason.

## 3. Expected outcomes on the 23 in-envelope points

Bands, adopted unchanged from the module scorecards so this exam cannot grade itself
leniently: **3.0x on a level row** (B2.1's `PASS_BAND_LEVEL`), **3.0x on acrylamide levels**
(B3's gating band). Fold error is reported for every in-envelope point regardless of band.

### 3.1 Hofmann 1998, 5 bundles / 10 points — expect **2/10 to 4/10**

Four of these ten points have already been scored, in B2.1, from the declaration's own row
inventory rather than from these files: `hofmann_ribose_pH3_MFT` 0.0589x FAIL,
`hofmann_ribose_pH3_FFT` 24.8x FAIL, `hofmann_ribose_pH7_MFT` 1.13x PASS,
`hofmann_ribose_pH7_FFT` 0.379x PASS. **Recorded here in advance, because it qualifies the
"first time any build wave opens these bundles" framing**: the FILES were never opened, but
the ribose pH-3 and pH-7 MEASUREMENTS were already scored under different row names. Both
sides were hold-out throughout, so nothing was fitted — but these four points are a re-score,
not a first exposure, and they must be labelled as such in the exam table.

Expectation: the two ribose pH-7 points pass again and the two ribose pH-3 points fail again
(the pH-3 rows are the module's declared weakest axis). The six genuinely new points are
glucose/cysteine at pH 3 and pH 7, and xylose/cysteine at pH 5. Glucose+cysteine reaches the
thiols only through the minor `r_glc_ha` / `r_glc_fur` fragmentation channels rather than the
pentose trunk, so **the glucose bundles are expected to under-predict hard, worse than the
ribose ones**. Xylose pH 5 is the closest to the fit condition (pH 5, 145 °C, 20 min) and is
the single most likely new pass.

### 3.2 Yiltirak 2026, 4 bundles / 8 points — expect **0/8 to 2/8**, under-predicting, worst at 100 °C

The genuinely new temperature axis, and the sharpest test in the exam. B2.1's constants are
`k_ref` at **145 °C** with a lumped formation Ea of 119.28 kJ/mol; these bundles sit at
100/110/120/130 °C, i.e. **15–45 °C BELOW the fit point**. Amendment 5 established that the
sulfur branch's real temperature behaviour is **non-Arrhenius** — apparent Ea rises from 7 to
98 kJ/mol across the legs — and B2.1's `kang_140C_MFT` row failed by 1500x extrapolating
*upward* over a 20 °C gap. A single lumped Ea extrapolating *downward* over a larger gap
should fail in the mirror direction.

Direction is pre-registered: **the model should OVER-predict at the low-temperature end** if
the true Ea is larger than the lumped 119.28 (cooling the model too little), and
**under-predict** if smaller. Kang's measured behaviour — the sulfur class decelerating while
apparent Ea rises — implies the true low-T Ea is *smaller* than the fitted lump, so the
pre-registered call is **under-prediction that worsens as temperature falls**, with the
100 °C / 4 h bundle the worst point in the family. If instead the errors are flat across the
ladder, that falsifies the lumped-Ea diagnosis and is the more interesting result.

### 3.3 Acrylamide, 5 points (Chang ×2 + Chang-water + Ye + Lin) — expect **0/5 to 1/5**

B3 scored **0/4** on its own gating rows and its one genuinely open row,
`knol2010_molar_yield_on_asparagine`, came in **3.84x LOW** against a 3.0x band. These five
points are all at **180 °C**, against a module whose fitted `k_ref` sits at 160 °C and whose
`Ea_int1_mel` hit the **top of its search bound at 260.0 kJ/mol** — a saturated bound is a
declared underfit, and extrapolating a saturated barrier 20 °C upward amplifies it.

Pre-registered direction: **under-prediction**, consistent with Knol 2010. The Lin 2022 point
is fructose-fed and therefore the weakest of the five — the module carries **no
fructose-specific chemistry at all** (B3's own words), reaching acrylamide only via the
trunk's `r_fru_glc` isomerisation, and B3's `dv1_fructose_over_glucose_acrylamide` diagnostic
already read 0.13x. Expect Lin to be the worst acrylamide point.

The three Chang/Ye glucose points differ only in time (10 vs 30 min) and loading (27.75/33.3
vs 200/200 mM). The 200 mM Ye bundle is a **7x higher loading** than Chang; because
initiation here is genuinely bimolecular, the predicted molar yield should RISE with loading.
That is the same axis `knol2010_molar_yield_on_asparagine` tests, and it is the one place in
this exam where the core's architecture makes a prediction the old lane's fixed-yield
fraction could not make at all. **Whether Ye/Chang moves in the right direction is a more
informative result than whether either passes.**

### 3.4 Summary expectation, stated as a number so it can be wrong

**Pre-registered overall: 2 to 7 of the 23 in-envelope points inside band, most likely 4.**
Median fold error pre-registered as **10x–100x**, i.e. NOT better than the old lane's
maillard_path median of 10.8638x, and plausibly worse. The core is not expected to beat the
old lane on median accuracy in this exam. What it is expected to do differently is (a)
refuse 17 points instead of guessing them, and (b) fail with named, localised causes.

## 4. What would falsify the cutover

Stated in advance so it cannot be rationalised afterward:

- If the core scores **0/23**, the cutover is wiring that shipped a worse predictor, and the
  README must say so in those words.
- If the core beats the old lane on the maillard_path median, that is a **positive** result
  and must be reported without inflating it — the old lane's 10.8638x median is itself over
  a different point set (it scored HMF and DMHF points the core declines).
- **The two medians are not directly comparable**, and the exam report must print the
  paired-subset median (old lane restricted to the 23 points the core answers) alongside the
  headline. A comparison over different denominators is the error this repository has already
  made twice.
- If any FAST-lane absolute ppb survives in a user-facing surface after the cutover, the
  cutover is incomplete regardless of the exam score.

## 5. Standing rules for this wave

1. No parameter, bound, prior or initialisation changes. The exam is run **once**.
2. If the exam reveals a WIRING bug, the wiring is fixed, the fix is logged in the exam
   report with a before/after, and the exam is re-run. A wiring fix is not a tuning.
3. Out-of-envelope points stay out. Widening the envelope after seeing a score is the
   prohibited move.
4. Every number in the exam report carries its lane, its band, and whether it is a re-score
   or a first exposure.
