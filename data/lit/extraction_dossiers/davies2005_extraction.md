# Davies 2005 — EXTRACTION (glucose 0.2 M + glycine 0.2 M in sodium acetate 0.1-0.5 M, initial pH 5.5, 50 °C, A420 to 2.0 over up to ~500 h: the late-stage browning rate rises ~6× from 0.1 to 0.5 M acetate, close to first order in acetate, no saturation; glucose + 0.5 M acetate browns WITHOUT glycine, glucose + 0.5 M citrate does not)
DOI: 10.1021/bk-2005-0905.ch012 (printed in the page margin; Crossref returns "Effect of Acetate Concentration on Maillard Browning in Glucose-Glycine Model Systems", Davies, ACS Symposium Series pp. 146-156, 2005-07-28; checked 2026-11-01)

**Read 2026-11-01 from `data/articles/davies2005.pdf`** (11 pp., ACS Symp. Ser. 905, Process and Reaction Flavors,
Weerasinghe & Sucan eds., chapter 12, pp. 146-156; an ABBYY FineReader scan). All pages read from the page images;
Tables 1-4 read from the images (the OCR prints "±" as "∀"/"V"). **Figs. 1-3 read from the images at screen
resolution** (read error ±0.05 A420, ±10 h). Table 3's powers of ten are cut off in the scan itself [!].

**Legend.** **[V]** verbatim · **[M]** printed number · **[F]** read off a figure · **[Q]** qualitative claim · **[C]**
cited · **[P]** proposed, not shown · **[Z]** derived by me · **[!]** flag.

## 0. Identity

| field | value |
|---|---|
| Title | "Effect of Acetate Concentration on Maillard Browning in Glucose-Glycine Model Systems" |
| Author | Cathy Davies — Department of Animal and Food Sciences, University of Delaware, Newark, DE 19716 (sole author) |
| Venue | ACS Symposium Series 905 (2005), chapter 12, pp. 146-156; published 2005-07-28 |
| Paper type | browning kinetics (A420 vs time), empirical regressions (linear late-stage slope; A = a t² + b t) |
| Lab | see sec. 4 — Delaware, with a Leeds design lineage |

## 1. What this paper answers for the model

The question (E5; Leeds as a development lab in `docs/reference/contrast_benchmark/split.yml`): the carboxylate
catalysis law; is this the Leeds lab?

| what the model needs | what Davies 2005 gives | tag |
|---|---|---|
| rate vs acetate concentration at fixed pH | **k0 (late-stage A420 slope) = 1.24 / 2.26 / 2.73 / 6.44 / 7.72 × 10⁻³ h⁻¹ at 0.1 / 0.2 / 0.3 / 0.4 / 0.5 M acetate** (±0.07-0.54; Table 1); ×6.2 from 0.1 to 0.5 M | [M] |
| order in acetate | log-log slope **1.15** on k0, **1.30** on the t² coefficient a; k0 / [Ac] = 12.4 / 11.3 / 9.1 / 16.1 / 15.4 × 10⁻³ h⁻¹ M⁻¹ — **no saturation to 0.5 M**, if anything slightly super-linear | [Z] |
| zero-buffer control | glucose + glycine with no acetate: ΔA420 **0.03 in > 300 h**, "essentially no reaction" [Q]; but its pH was NOT adjusted: **initial 5.9, final 4.9** vs 5.5 → 5.3 in every acetate pot [M] | [M] / [!] |
| pH during the run | buffered pots 5.5 → 5.3 [M]; acetate base fraction at pH 5.5 ≈ 0.85 (pKa 4.76) [Z] | [M] / [Z] |
| acetate vs citrate | amine-free (0.2 M glucose, no glycine): **0.5 M acetate browns to A420 ≈ 0.6 (plateau ~350-500 h); 0.5 M citrate ≈ 0 to ~280 h** (Fig. 3) | [F] |
| acetate vs phosphate | **not tested** | — |
| temperature | **one temperature, 50 °C** | [M] |
| sugar and amine orders | doubling glucose OR glycine at 0.3 or 0.5 M acetate raises a ~3× (Table 4) [M/Q]; "regardless of whether" glucose or glycine was increased [Q] | [M] |
| intermediates | none measured | — |

**Answers.**

1. **The carboxylate law here is k0 ≈ k_Ac [acetate]^1.1-1.3, with a near-zero intercept** (confidence in the form
   ~0.6). Linear fit k0 = 17.1 × 10⁻³ [Ac] − 1.1 × 10⁻³ h⁻¹ [Z]. There is no saturation between 0.1 and 0.5 M. The
   no-acetate pot is not at the same pH (it starts higher, 5.9, and ends lower, 4.9), so the "×∞" contrast against
   zero buffer is not clean; the 0.1-0.5 M ladder is (pH 5.5 → 5.3 in all).
2. **The author's reading** — acetate takes part directly, "not caused by improved pH control" [V] — is [P]: no
   mechanism is shown. The amine-free browning with acetate (and not with citrate) at 50 °C is the strongest datum
   for an acetate-specific sugar-degradation path [F].
3. **Same lab as Leeds? Not established by the paper** (sec. 4). On the conservative merge rules it should probably
   count with `leeds_wedzicha` (owner decision).

## 2. Conditions

- Glucose 0.0-0.4 M, glycine 0.0-0.4 M (standard 0.2 + 0.2 M); acetate buffer **0-0.5 M, pH 5.5, adjusted with glacial
  acetic acid** (sodium salt implied, counter-ion not stated) [M].
- **50 °C**, incubated, sampled at timed intervals until A420 reached 2.0; up to ~500 h [M/F]. Vessel, headspace and
  O2: **not stated** [!]. No chelator, no azide stated.
- pH: Corning 443i meter, at room temperature presumably (not stated) [M]. Water activity: dilute aqueous.
- Analysis: A420 on a Beckman DU60 against a water blank; triplicates ("All reaction mixtures were prepared in
  triplicate") [M]. Regressions: Excel (linear), NLREG 5.2 (A = a t² + b t).

## 3. Quantitative results

**Table 1** (0.2 M Glc + 0.2 M Gly, pH 5.5, 50 °C) [M]:

| acetate M | 10³ k0, h⁻¹ | intercept c | r² |
|---:|---:|---:|---:|
| 0.1 | 1.24 ± 0.07 | −0.162 ± 0.01 | 0.9860 |
| 0.2 | 2.26 ± 0.2 | −0.241 ± 0.03 | 0.9744 |
| 0.3 | 2.73 ± 0.2 | −0.300 ± 0.03 | 0.9787 |
| 0.4 | 6.44 ± 0.5 | −0.709 ± 0.04 | 0.9815 |
| 0.5 | 7.72 ± 0.5 | −0.721 ± 0.06 | 0.9798 |

[Z] lag (x-intercept −c / k0): 131 / 107 / 110 / 110 / 93 h — the induction period shortens only mildly with acetate.

**Table 2** (A = a t² + b t) [M]: 10⁶ a = 2.84 / 5.41 / 6.88 / 15.7 / 24.7; 10⁴ b = −2.0 / −2.0 / −3.0 / −5.0 / −10.0
(0.1 → 0.5 M); r² 0.970-0.999. The abstract restates these.

**Table 3** [M, exponents illegible]: a = m[Ac] + n[Ac]²; b = p e^(−q[Ac]); m = 8.99 ± 1.98, n = 3.23 ± 8.35, p = 5.64
± 1.78 (all "× 10^?", the exponent cut off in the scan), q = −5.71 ± 0.69. **[!]** n's error exceeds its value; with q
negative and b negative, p must be negative for the fit to reproduce Table 2, but it prints positive.

**Table 4** (10⁴ × a and b; glucose / glycine / acetate) [M]:

| Glc M | Gly M | b, 0.5 M Ac | b, 0.3 | b, 0.1 | a, 0.5 | a, 0.3 | a, 0.1 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.4 | −32.1 | −13.6 | 6.0 | 0.85 | 0.27 | 0.04 |
| 0.2 | 0.2 | −10.0 | −2.8 | −1.5 | 0.24 | 0.07 | 0.03 |
| 0.4 | 0.2 | −49.3 | −8.1 | 2.3 | 1.04 | 0.22 | 0.05 |

(The 0.2 / 0.2 / 0.5 M row matches Table 2: a = 24.7 × 10⁻⁶, b = −10 × 10⁻⁴ — consistent.)

**Fig. 1** [F]: 0.5 M acetate, A420 ≈ 0 until ~80 h, then 0.12 (125 h), 0.65 (200 h), 1.37 (270 h); fitted line y =
0.0077x − 0.7212, R² 0.9801 (printed on the figure) [M].

**Fig. 3** [F] (0.2 M glucose, pH 5.5, 50 °C): acetate 0.5 M + glycine 0.2 M → 1.37 at ~270 h; acetate 0.5 M, no
glycine → 0.25 (~150 h), 0.40 (250 h), 0.55-0.62 (380-500 h, plateau); citrate 0.5 M, no glycine → ≤ 0.02 to ~280 h.

## 4. Mechanism evidence and its kind

- **Within-study factor** (kind 6): acetate concentration ladder at a held start pH; glucose / glycine doublings; a
  buffer-identity contrast (acetate vs citrate) without amine. Fit only on the late-stage slope; no intermediate,
  no product identity.
- The proposed 3-step scheme A + S → I1 → I2 → B (+A), B = k'[A]²[S]t² is cited from Davies, Wedzicha & Gillard 1997
  [C]; the acetate terms are empirical (Table 3) [P].
- **Lab.** Single author at the University of Delaware. Reference 1 is "Davies, C. G. A.; Wedzicha, B. L.; Gillard,
  C. Food Chemistry 1997" (the Leeds glucose-glycine kinetic model) and reference 13 is Swales & Wedzicha 1995
  (acetate and sulfite loss, Leeds); the design (glucose-glycine, acetate buffer pH 5.5, 50-55 °C, A420) is the
  Leeds design. **Whether "Cathy Davies" is C. G. A. Davies of the Leeds paper is not stated** [!]; I think it likely
  (~0.8) [Z]. By `lab_ids.yml` rule 3 a junior author who moved does not carry the lab, which would make this a new
  (TEST) laboratory; by rule 2's conservatism ("over-merging can only move cells from TEST to DEVELOPMENT") and the
  reused design, merging with `leeds_wedzicha` (a DEVELOPMENT lab) is the safe choice. **Owner decision.**

## 5. Bearing on open questions

- **E5 (no buffer-species catalysis): SUPPORTS** for carboxylate catalysis at 50 °C (×6.2 over 0.1-0.5 M acetate at
  held pH). **CONTRADICTS "saturating"** in the proposed fix: order ~1.1-1.3 in acetate, no plateau to 0.5 M
  (agrees with Potman 1989's order 1.34 in phosphate, `potman1989_extraction.md`). Leaves "fading at pH 8" untested
  (one pH).
- **A transfer warning for the Leeds constants** [Z]: the Leeds kinetic pots are in acetate (Wedzicha & Edwards 1991:
  50 mM acetate pH 5.5, `wedzicha1991_extraction.md`; Mundt & Wedzicha 2004: 0.2 M acetate pH 5.5,
  `mundt2004_extraction.md`). If acetate's effect is general (not browning-specific), Leeds rate constants carry a
  buffer term of roughly ×2 per doubling of acetate and should not be transferred to unbuffered or phosphate pots
  without a declared band. Which steps acetate acts on is not shown here.
- **CF-011: LEAVES UNCHANGED** (no phosphate, no cysteine). It weakens "phosphate-specific" readings in general:
  a carboxylate at 0.1-0.5 M catalyses Maillard browning comparably in order.
- **F14 amine-free degradation:** acetate (0.5 M, pH 5.5, 50 °C) degrades glucose to colour without amine; citrate
  does not — a buffer-identity effect at low temperature [F].
- Evidence map 1.2: one more strict within-study buffer ladder (acetate), one lab (Leeds lineage).

## 6. Disagreements found

1. **Citrate amine-free:** none at pH 5.5 / 50 °C here vs ×12 over no buffer at pH 3 / ~121 °C in Saunders & Jervis
   1966 (`saunders1966_extraction.md`, Table II). Causes: temperature (carboxylate catalysis or citrate decomposition
   needing ≥ 100 °C), pH, and response (A420 at 50 °C vs T490 at 121 °C).
2. **Acetate per mole vs other anions:** here acetate alone gives large effects at pH 5.5; Rizzi 2004 puts acetate at
   ~1/8 of phosphate's A420 effect (pH 7.8 start, unbuffered acetate, ~100 °C, `rizzi2004_extraction.md`); Saunders
   puts it at ~1/3 of phosphate / citrate per mole at pH 3. Causes: acetate's base fraction (0.85 at pH 5.5 vs ~1 at
   7.8 vs 0.02 at 3.0 [Z]), temperature, sugar.
3. **The amine-free acetate plateau** at A420 ≈ 0.6 (Fig. 3) with no glycine is unexplained: an exhausted minor
   reactant (an impurity in the acetate?) would give such a plateau; the paper does not discuss it [!].
4. **Table 3's fitted p has the wrong sign** to reproduce Table 2's negative b (sec. 3) [!].

## 7. Ledger-ready measurements

1. **Table 1, all five rows** (printed): k0 vs acetate at held pH 5.5, 50 °C — the cleanest carboxylate
   concentration law on disk; key k0 with its SD and the lag [Z] as a note.
2. **Fig. 3** (two readers): glucose + 0.5 M acetate vs + 0.5 M citrate, no amine — a buffer-identity contrast at
   50 °C.
3. Table 4 (printed): glucose and glycine doublings at three acetate levels — the reactant orders at fixed buffer.
4. Not worth keying: the no-acetate pot (unadjusted pH, no rate printed beyond "0.03 in > 300 h").
