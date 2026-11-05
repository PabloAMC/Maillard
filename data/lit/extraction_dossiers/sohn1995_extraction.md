# Sohn & Ho 1995 — EXTRACTION
DOI: 10.1021/jf00060a001 (the PDF prints the manuscript code "JF950087L"; Crossref 2026-11-05: "Ammonia Generation during Thermal Degradation of Amino Acids", Sohn, Ho, J. Agric. Food Chem. 43(12), 3001-3003, 1995. **Matches the brief and the PDF.**)
### Ammonia released by 19 amino acids heated ALONE in water (no sugar), initial pH 8 (NaOH, unbuffered), 2 h in a stainless-steel reactor, ammonia-selective electrode. **Cysteine: ~82 % of its α-amino N released as NH3 at 180 °C** (0.1 M-basis run and 0.05 M run agree); **temperature series (0.05 M, 2 h): ~14 % at 110 °C, ~76 % at 150 °C, 82 % at 180 °C** [F]. Glutamine 100 % (amide) at every T; glutamic acid ~1 %. **One time point per temperature, no cysteine or product assay, no pH after heating.**

**Read 2026-11-05 from `data/articles/sohn1995.pdf`** (3 pp., J. Agric. Food Chem. pp. 3001-3003; scan with a text
layer). Prose from `pdftotext -layout`; **Figs. 1 and 2 digitised from 250 dpi renders** (scratchpad
`img/sohn2f1.png`, `sohn3f2.png`). Read error ± 0.002 mol NH3 per 0.1 M (Fig. 1) and ± 0.0015 per 0.05 M (Fig. 2)
[F]. The only printed cysteine number is the text's 82 % at 180 °C.

**Legend.** **[V]** verbatim (quote < 15 words) · **[M]** number printed in this paper · **[F]** read off a figure ·
**[Q]** qualitative claim by the authors · **[C]** cited · **[P]** proposed, not shown · **[Z]** derived by me, never
printed · **[!]** flag.

**Laboratory.** Department of Food Science, Rutgers (Ho). **`ho_network` in `lab_ids.yml`** (the 1980s-90s Rutgers Ho
papers join through Ho himself), which is one of E10's three DEV laboratories.

## 1. What this paper answers for the model

The case (M1c sec. 4.3 and 5.2) uses it as a band source for R-sink's cysteine-only arm, and for CF-084 (M1b: the
engine's charge balance declares the cysteine amine destroyed, `ph_state.AMINE_FATE_BASIS`).

| model need | what Sohn 1995 gives | tag |
|---|---|---|
| cysteine-alone degradation in water vs T | NH3 yield after 2 h, initial pH 8: ~14 % (110 °C), ~76 % (150 °C), 82 % (180 °C) | [F] / [M] |
| fate of the cysteine amine | released as free ammonia, not destroyed (the majority at ≥ 150 °C) | [M] / [Q] |
| a rate constant | **none printed**; one time point per T | [Z] only |
| the sulfur fate (H2S, cystine), cysteine left | **not measured** | — |
| other amino acids, same design | Asn 0.125 / 0.1 M (≥ 100 % incl. amide), Asp ~60 %, Gln ~95-100 %, Arg ~14 %, others ≤ 10 % at 180 °C | [F] |

## 2. Conditions

| item | value |
|---|---|
| Fig. 1 runs | 2 g amino acid in 250 mL deionised distilled water (cysteine: 8 g/L = **66 mM**), results "calculated on the basis of 0.1 M"; initial pH adjusted to 8 with 0.1-1 N NaOH; **no buffer** |
| Fig. 2 runs | Asn, Asp, Gln, Glu, Cys at **0.05 M** each, pH 8, **110, 150, 180 °C**, 2 h |
| Vessel | 0.3 L stainless-steel reactor (Hoke), headspace ~50 mL if the 250 mL charge was transferred whole (my reading; atmosphere not stated, so air) |
| Heating | conventional drying oven preset to temperature (heat-up time not stated); 2 h; ice-bath quench |
| Assay | Orion 95-12 gas-sensing ammonia electrode, 5 mL ionic-strength adjustor per sample, NH4Cl standards 10⁻¹-10⁻⁴ M on a semilog calibration per batch |
| Replicates | "a minimum of triplicates"; error bars on the figures (SD presumably; not stated) |
| pH after heating | **not reported** |

## 3. Figures digitised

### Fig. 1 (p. 3002). mol NH3 per 0.1 M amino acid, 180 °C, 2 h, pH 8 [F]

| Ala | Arg | Asn | Asp | **Cys** | Gln | Glu | Gly | His | Ile | Leu | Lys | Met | Pro | Ser | Thr | Trp | Tyr | Val |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ~0 | 0.014 | 0.125 | 0.060 | **0.082** | 0.095 | ~0 | ~0 | 0.006 | 0.004 | 0.003 | 0.005 | 0.004 | ~0 | 0.006 | 0.004 | 0.010 | 0.008 | ~0 |

Text percentages [M]: Arg "about 28 %" (attributed to the guanidino group), His 6.2 %, Trp ~10 %, Ser / Thr 5-6 %,
non-polar < 5 %, Gln 100 %, Glu 1.3 %. The Arg bar (0.014 / 0.1 M) and the printed 28 % do not match on a 0.1 M basis
(Flag 3).

### Fig. 2 (p. 3003). mol NH3 per 0.05 M amino acid, 2 h, pH 8 [F]

| amino acid | 110 °C | 150 °C | 180 °C |
|---|---|---|---|
| **cysteine** | **0.0070** (14 %) | **0.0381** (76 %) | **0.0410** (82 %; text: 82 % [M]) |
| asparagine | 0.002 | 0.029 | 0.063 |
| glutamine | 0.048 | 0.051 | 0.050 |
| aspartic acid | ~0 | 0.003 | 0.030 |
| glutamic acid | ~0 | ~0.001 | ~0 |

Cysteine error bars: visible only at 180 °C (± ~0.001).

## 4. Arithmetic (mine) [Z]

**1. First-order equivalents of NH3 release (cysteine, 0.05 M, pH 8 initial).** k = −ln(1 − f)/120 min:
110 °C **1.26 × 10⁻³ min⁻¹**; 150 °C 1.20 × 10⁻² min⁻¹; 180 °C ≥ 1.4 × 10⁻² (yield plateau, not usable). Apparent Ea
from 110 and 150 °C: **~76 kJ/mol** (two points, the 150 °C one at 76 % conversion; heat-up time unknown, which
biases the 110 °C point low). Extrapolated: ~3.4 × 10⁻⁴ (90 °C), ~6.6 × 10⁻⁴ (100 °C), ~2.3 × 10⁻³ min⁻¹ (120 °C).

**2. What NH3 measures.** NH3 is a LOWER bound on deaminated cysteine (ammonia partitions to the headspace at pH 8
and can re-add to carbonyl products) and it is not a cysteine-loss rate: cysteine lost to cystine or to adducts
without deamination is invisible. Equally, one cysteine releases at most one NH3, so the 82 % at 180 °C means ≥ 82 %
of the cysteine was degraded there.

**3. The plateau.** 150 → 180 °C adds only 6 points (76 → 82 %); 18 % of the N is in a form that does not yield NH3
(amine-retaining products, or NH3 lost from the liquid).

## 5. What it can and cannot test in the model

- **Can:** give a cysteine-ALONE (no sugar) deamination yield vs temperature at initial pH 8 in water, from an
  E10 DEV laboratory; give a first-order-equivalent order of magnitude for a cysteine-only sink at 110 °C
  (≥ 1.3 × 10⁻³ min⁻¹, on NH3) with a crude Ea; show that the cysteine amine is largely RELEASED AS AMMONIA, a base
  (pKa 9.25), not destroyed: relevant to CF-084's amine-fate basis for a cysteine-only degradation.
- **Cannot:** give a rate (one time per T, unknown heat-up), the pH trajectory (not measured), the sulfur fate, the
  cysteine left, or anything at the Maillard pots' acid pH (4-6.5). Stainless steel and air may catalyse thiol
  oxidation (not tested).

## 6. Conflicts with existing records

| record | what it says | what Sohn says | verdict |
|---|---|---|---|
| `zhou2022_extraction.md` (sec. Arithmetic: Cys 20 mM alone, 120 °C, label pH 7: ~2.5 × 10⁻³ min⁻¹ from Cys loss) | Cys loss rate | NH3-release equivalent ~2.3 × 10⁻³ min⁻¹ at 120 °C by extrapolation | **consistent in order of magnitude** (same lab network; different observable and pH 8 vs 7) |
| `zhao2017b_extraction.md` / `cao2017_extraction.md` (all cysteine-destroying channels ≤ 2-4 × 10⁻⁴ min⁻¹, 90 °C, pH 4.1-4.7, with xylose) | ceiling | ~3.4 × 10⁻⁴ min⁻¹ at 90 °C by extrapolation, pH 8, no sugar | no conflict: pH 8 vs 4.3 and a 2-point Ea; at the ceiling's edge |
| `mizutani1960_extraction.md` (Hokkaido; free cysteine, 0.1 M, water, 120 °C oil bath, CO2 sweep: NH3 0.75 % of N at 3 h, 1.34 % at 6 h; H2S 0.96 % of S at 6 h) | ~4 × 10⁻⁵ min⁻¹ NH3-equivalent at 120 °C [Z] | ~1.3 × 10⁻³ min⁻¹ at 110 °C | **conflict of ~30-60x in deamination rate** (Sohn's 110 °C value is already 30x Mizutani's 120 °C one): pH 8 (NaOH) in a sealed stainless-steel reactor vs cysteine's own pH under a CO2 sweep (open, acidic). It argues the cysteine-only sink is strongly base-favoured, or vessel-catalysed; neither is separated by either paper. `obata1960_extraction.md` (same lab, 100 mL) measures H2S only. (Mizutani's F04#26 is E10's TEST contrast; this is a measured-vs-measured comparison, no engine.) |
| CF-084 (the engine declares the cysteine amine destroyed in `ph_state`) | amine destroyed | amine released as NH3 (≥ 76 % at ≥ 150 °C, 2 h; 14 % at 110 °C) | **bears on CF-084** for cysteine-alone degradation only; the Maillard pots' amine fate (Amadori, TTCA, pyrazines) is not measured here |

## 7. Flags

1. **[!] One time point per temperature; oven heating with no heat-up time** — not a kinetic study.
2. **[!] pH 8 is the initial, unbuffered pH**; NH3 release would raise it, acid products lower it; neither measured.
3. **[!] Fig. 1 basis.** 2 g in 250 mL is 66 mM for cysteine and 30-150 mM across the set, but results are "on the
   basis of 0.1 M"; the conversion is not shown. The Arg text (28 % from the guanidino group) is not reconcilable with
   its bar on that basis. Fig. 2 (0.05 M, consistent with the text's 82 %) is the cysteine number to use.
4. Ammonia electrode reads free NH3 + NH4⁺ after the ISA makes the sample alkaline; volatile losses before reading
   are not quantified.

## M1c role (frozen bd641ae8)

**Role fixed by the case (sec. 4.3 and 5.2):** band source for R-sink's cysteine-only arm (by amendment), and for
CF-084 (M1b).

**Band it can supply:**
- **Coordinate: R-sink's cysteine-only sink** (case band 0-2 × 10⁻⁴ min⁻¹ at 90 °C primary; at 120 °C, Zhou 2022's
  ~2.5 × 10⁻³ as the upper corner). Sohn gives, for cysteine alone in water at initial pH 8: a lower bound of
  **≥ 1.3 × 10⁻³ min⁻¹ at 110 °C** (NH3-release equivalent), apparent Ea ~76 kJ/mol (110-150 °C), extrapolating to
  ~3.4 × 10⁻⁴ (90 °C) and ~2.3 × 10⁻³ min⁻¹ (120 °C). The 120 °C value independently supports Zhou 2022's upper corner
  (same lab network); the 90 °C value sits just above the primary band's top, at pH 8 rather than the pots' 4-6.5.
  Mizutani 1960's cysteine at its own (acid) pH deaminates ~30-60x slower at 120 °C, so at the pots' pH the sink is
  probably well below Sohn's pH-8 numbers: Sohn is an upper-corner source, not a centre.
- **Transfer:** cysteine itself, but no sugar, initial pH 8, stainless steel, NH3 not Cys; magnitude good to a factor
  ~2-3 at best.
- **CF-084 (M1b):** in cysteine-alone degradation the amine leaves largely as ammonia (a base), 14 % (110 °C) to 82 %
  (180 °C) in 2 h; this argues against declaring that amine destroyed for the cysteine-only sink's share.
- **Not decisive.** A cysteine-only sink cannot make the 1 : 1 Cys : Xyl deficit (case sec. 2.9 and 4.2), so by the
  case's own logic nothing here can be decisive for (c).
