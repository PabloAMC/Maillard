# Wave B48 pre-registration — the 2,4-decadienal breakage step on the lipid lane

*Written 2026-10-27 before any B48 number existed. Code: `src/kinetic_core/parameters_lipid.py`
(`DECADIENAL_BREAKAGE`), `src/kinetic_core/lipid.py` (`integrate_lipid`); ship rule:
`scripts/generators/generate_kinetic_core_b48_ship_rule.py`.*

## 1. Why

The lipid lane declares `DECADIENAL_TO_HEXANAL_GAP`: "No rate for it exists in the corpus." That stopped
being true on 2026-09-11, when wave B42 read Zamora, Navarro, Aguilar & Hidalgo 2015 (Food Chem 174:89,
`zamora2015b_extraction.md`) and named the step as a candidate; the gap text was never updated. On
2026-10-27 Zamora, Gallardo & Hidalgo 2007 (`zamora2007_extraction.md`) gave an independent number from
another experiment of the same laboratory. B48 builds the step from the printed numbers. **Nothing is
fitted.**

## 2. The step

2,4-Decadienal breaks thermally, first order in itself, through two measured channels:

| channel | printed yield, 1 h at 200 °C | Ea (initial rates, 120/160/200 °C) | as a first-order constant at 200 °C |
| --- | --- | --- | --- |
| → hexanal (+ a C4 fragment) | 11.5 % (slope 0.1154, r > 0.998) | 21.3 kJ/mol | −ln(1 − 0.1154)/60 = 2.045e-3 /min |
| → 2-octenal (+ a C2 fragment) | 0.8 % (slope 0.00821) | 29.6 kJ/mol | −ln(1 − 0.00821)/60 = 1.374e-4 /min |

`k(T) = k(200 °C) exp(−Ea/R (1/T − 1/473.15 K))`. The lane solves each segment analytically, and B48
keeps that: decadienal made from the hydroperoxides during a segment and broken within it is the
exact two-step first-order solution, not a lag. Hexanal is credited; 2-octenal is not a lane species,
so its carbon and both fragments go to `LIPID_FRAG_C`, and the carbon balance still closes as an
equality.

**Declared limits.** The source medium is tetrahydrofuran–buffer, pH 8 (the yield peaks there), closed
tube, 4 µmol in 0.5 mL: not a food. Other losses of decadienal (the remainder of its disappearance,
figure-only in the source) are not modelled, so the decadienal level is an upper bound on what survives.
The step adds no pH term. **Window:** the barriers come from 120/160/200 °C, so the step is OFF below
100 °C (one 20 °C rung under the measured range) rather than extrapolated: a 21 kJ/mol barrier carried to
storage would break most of the decadienal in weeks at room temperature, which no source says.

**Does it double-count B6?** No. B6's fit column is Frankel & Gardner 1989, hydroperoxides thermolysed
in a GC injector at 180 °C: seconds, against a half-life of tens of minutes. The breakage was not running
during the measurements B6 fitted, so B6's split stands and is not refitted.

## 3. Expectations, stated in advance

- At 140 °C for 30 min the hexanal channel converts ~2.6 % of the decadienal present; at the UHT and
  40 °C pots on the panel it is negligible. Scored lipid rows should move by well under 0.1 dex.
- Zamora 2007 (pH 3, phenylalanine present, 180 °C, 1 h) printed 9–14 % of decadienal to hexanal. The
  step, built from Zamora 2015b alone, predicts ~9 %: inside the range at its low end.

## 4. Ship rule

- **T1 (hold-out):** pure 2,4-decadienal held 60 min at 180 °C: the predicted share converted to hexanal
  lies within Zamora 2007's 9–14 %, widened by 1.5x on each side for the medium (pH 3 vs 8, amino acid
  present): [6 %, 21 %].
- **T2:** the scored lipid rows (trust loop and external matrix): within-3x count not lower than without
  the step; median |log10 error| not worse by more than 0.05 dex.
- **T3:** the directional panel, independent claims agreeing: not fewer than without the step.
- **T4 (recorded):** every scored row that moves by more than 0.01 dex.

**SHIP if T1, T2 and T3 hold**: the step is enabled and `DECADIENAL_TO_HEXANAL_GAP` is rewritten as
closed. Otherwise the step stays in the code disabled, with the failing test named.

## 5. What B48 does not build: the lipid-driven Strecker step

Reviewed the same day and **not built**. Zamora 2007, 2008, 2013 and 2015 and Hidalgo 2004 measure
Strecker degradation by lipid carbonyls for **phenylalanine only**. The Strecker aldehydes that matter
for meat aroma are methional, 2- and 3-methylbutanal (methionine, isoleucine, leucine); the model
refuses all three today (B22 found the identity assumption between amino acids refuted by two
laboratories). A lipid-driven step would therefore need either a target the model does not carry
(phenylacetaldehyde) or a transfer from phenylalanine to amino acids where B22 showed transfer fails.
What would unlock it: one paper measuring a lipid carbonyl (2,4-decadienal or a 4,5-epoxyalkenal)
with leucine or methionine at cooking temperature, with yields.

## 6. Verdict

*(Written after the run, below this line.)*

**SHIP** (`kinetic_core_b48_ship_rule.md`, 2026-10-27). T1: pure 2,4-decadienal, 60 min at 180 °C,
9.2 % to hexanal, inside Zamora 2007's printed 9–14 % (the step was built from Zamora 2015b alone).
T2: the eight scored lipid rows 4/8 within 3x with and without the step, median |log10 error| 0.448
both ways. T3: directional 25/43 both ways. T4: one row moved more than 0.01 dex, and it improved:
Bi 2020 roasted-pea hexanal, 3.66x → 3.39x — the only scored pot hot and long enough for the step to
act. The expectation in §3 held: the step is real and small on the present panel. It ships enabled;
`DECADIENAL_TO_HEXANAL_GAP` now reads as closed.
