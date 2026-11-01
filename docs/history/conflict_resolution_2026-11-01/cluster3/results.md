# Cluster 3 results (2026-11-01): in-situ pH and speciation

Scripts (all run through `./scripts/docker_maillard.sh run`, B51 as shipped, no engine switch changed, no refit):
`cf065_speciation.py`, `cf054_endpoint.py`, `ph_controls.py` (+ `ph_controls.log`), `ph_controls_inventory.py`,
`ph_arp_pka_sensitivity.py`, `buffer_capacity.py`, `ledger_scan.py` (+ `ledger_scan.log`). Prereg: `prereg.md`.

## 1. Controls of the pH model (run first; they decide the "systematic fix")

| control | data | shipped engine | verdict |
|---|---|---|---|
| positive, out of sample | Hou 2026 GX-ARP 5 mM, water, 7.0 -> 5.05 (20 min) -> 3.55 (120 min), 120 C | 6.35 at 20 min (+1.30), 4.26 at 120 (+0.71); mean +0.58, RMS 0.69. Sub-model 1 arm-B drift: RMS 1.09 | FAIL (tolerance 0.3) |
| same, ARP amine pKa set to 7.6 / 8.18 (measured DFG analogue) | as above | RMS 0.45 / 0.84; early too acid, late plateau: the engine's acid appears promptly with ARP loss (CBX carried out), the measured acid keeps coming after the ARP is gone | no single constant fits the SHAPE |
| within-lab | Zhou 2023 ARP 20 mM alone, labels 6 / 7 / 8 -> 3.98 / 4.11 / 4.95 | 3.19 / 4.55 / 7.09 (pKa 7.07) | label sensitivity ~4x too large |
| negative | Meynier 1995: Cys 50 + ribose 35 mM, 0.2 M pyrophosphate, 140 C 1 h, measured drift <= 0.2 | drift -0.79 / -0.65 / -0.39 at labels 4.5 / 5.5 / 6.5 (arm B: -1.1 / -1.25 / -0.69) | FAIL |
| diagnosis | end-state inventory | Meynier: CBX 37.6 mM (cysteine carboxyl carried out, amine declared destroyed, `AMINE_FATE_BASIS`), ACID x yield ~0 | the named failure direction of AMINE_FATE_BASIS ("too acidic") is observed in a buffered cysteine pot |
| buffer capacity (finite difference on the balance) | 0.5 M phosphate pH 5: beta 8.8 (25 C) / 12.2 (145 C) mmol/L/pH; 0.2 M pyrophosphate pH 4.5: 5.4 / 10.9; 0.2 M phosphate pH 5.5: 9.1 / 5.7 | -- | "buffered at pH 5" is a weak buffer; the triage's negative control premise ("0.5 M pots do not move") is wrong at pH 5 |
| in-situ label shift, Hofmann 0.5 M phosphate | labels 3 / 5 / 7 at 145 C | 3.71 / 5.27 / 7.24 at t = 0; the pH-5 pot then drifts to 4.47 in situ (cooled 3.93) in 20 min | already in the engine (B51) |

**Consequence:** scoring every unbuffered pot at the ENGINE-PREDICTED in-situ pH would add errors of 0.5-1.3 pH units
(larger than most label-to-in-situ shifts) and moves buffered cysteine pots the wrong way. A systematic fix must
use MEASURED trajectories as an operator, plus a declared band where none exists (sec. 6).

## 2. CF-065 (Hofmann T3 / T4 vs Zhai 2023 JAFC Fig. 6)

In-situ states (central): Hofmann pH 5.34 at 145 C (I 0.52 M), HS- 9.6 % of dissolved sulfide, dissolved fraction
0.15-0.24 of the charged 20 mM (Vg/Vl 3, Kaw 1.2-2.2); Zhai pH 6.33 at 120 C (the 20 mM NH4+ buffers it down from
the 7.0 label; NH4+ pKa 7.0 at 120 C), HS- 47 %, dissolved 0.34-0.70 (Vg/Vl 0.5-3). HS- ratio 4.9x, not the
triage's "~50x" (label-only ratio at temperature: 17x).

Implied second-order constant gap, log10(Zhai / Hofmann) at 120 C, median [5-95 %] over the declared bands
(Ea 40-100, Kaw clamp vs van 't Hoff, Vg/Vl, come-up, precursor means):

| correction | FFT | MFT |
|---|---|---|
| none (charged sulfide, no T) | +1.94 [1.61, 2.27] | +3.08 [2.61, 3.54] |
| + temperature (145 -> 120 C) | +2.49 [2.10, 2.90] | +3.63 [3.11, 4.15] |
| + headspace (dissolved sulfide) | +2.07 [1.61, 2.52] | +3.21 [2.64, 3.77] |
| + HS- speciation at in-situ pH | **+1.38 [0.94, 1.83]** | **+2.52 [1.96, 3.07]** |

P(gap <= 0.5 dex) = 0.00 for both. Negative controls pass (Hofmann T7 1 vs 2 mmol: operator 2.00x, measured 2.22x;
matched pair shift 0.000). Condition corrections remove ~1.1 dex of the temperature-corrected gap; 1.4 (FFT) and
2.5 (MFT) dex remain. The residual differs by 1.1 dex between the two thiols, so one common Zhai scale factor (R3)
cannot be the whole residual. R2 (NH3 co-reactant) predicts the wrong sign for yields on the only on-disk amine +
sulfide contrast (Zheng 1997: alanine + Na2S 39.6 vs Na2S alone 64.2 mg/mol sulfur volatiles, pH-confounded), but
Zhai's furfural is lost ~7x faster than FFT forms (fingerprint of a large non-thiol furfural sink, consistent with
NH3). Post hoc (NOT pre-registered): thiol loss during Hofmann's 145 C hold (MFT less stable than FFT) could account
for part of the MFT-FFT difference; Grosch 1992's in-pot MFT loss ~0.03-0.05 /min at 100 C bounds it only with an
unmeasured Ea.

**Verdict: EXPLAINED IN PART (COND ~1.1 dex), residual OPEN.** P(residual mainly OP, Zhai's RF = 1 scale) 0.45;
product loss in Hofmann's hold 0.25; NH3 chemistry 0.2; other 0.1. Decider: one fed pot, furfural or norfuraneol +
Na2S vs + (NH4)2S at matched in-situ pH, stated fill, SIDA after cooling.

**Model consequence:** the B2.3 fed rows `fed_furfural_h2s` / norfuraneol + H2S (`generate_kinetic_core_b2_3_fit.py`)
charge 20 mM H2S as dissolved; at Hofmann's Vg/Vl = 3 and 145 C only 15-24 % is dissolved, so `k_fur_fft_hs` and
`k_nf_mft` are biased low by 0.6-0.8 dex relative to small-headspace pots (Hofmann's own 100 mL cysteine pots:
dissolved ~34 %). Proposed engine change: an H2S vessel-partition operator (Kaw(T) from henry_constants.yml, Vg/Vl
per pot, declared band when unstated). Zhai Fig. 6 stays shape-only.

## 3. CF-054 (Chen 2000 vs Chen & Ho 2002, 180 C, final pH)

SID from the label: +64.1 mmol/L NaOH at pH 8.5, ~0 at pH 5. Acid needed for each endpoint over bands NH4+ 0-60,
H2S 0-40, residual Cys 0-30 mM:
- start 8.5: 6.05 vs 4.96 needs an acid ratio of 1.23 [1.20, 1.31]; dpH/dln(acid) = -4.35 (endpoint on the
  equivalence slope).
- start 5.0: 4.01 vs 3.86 needs 1.18 [1.17, 1.19] wherever any NH4+ is released (the NH4+ = 0 corners, which need
  < 1 mM acid in total, are unphysical); dpH/dln(acid) = -0.92.
- One ~20 % difference in total acid between the two runs reproduces BOTH pairs. With the 10 % replicate CV (Davidek
  2002), the sd of a two-run difference is 0.61 pH (start 8.5) and 0.13 (start 5); observed 1.09 (1.8 sd) and 0.15
  (1.2 sd).

**Verdict: EXPLAINED (E1, run-to-run acid spread amplified ~5x by the alkaline start), P 0.75.** Consequence: none
for the engine; any use of an alkaline-start unbuffered endpoint carries ~0.6 pH units of irreducible error.

## 4. CF-069 (r_aa_frag)

- The trunk lane carries NO pH state (`ph_state.UNTRACKED_TITRATABLE["Gly"]`), while the fitted Martins 2005 pots
  (0.1 M phosphate, beta ~40-47 mmol/L/pH) drift 6.76 -> 5.51 in 240 min at 100 C (martins2005b ladder, measured).
- Sub-model 1 (inline pH, on disk): arm B (k_aa_frag OFF) fits Martins 2005 AA at RMS 0.16 vs arm A (free) 0.14, and
  predicts the held-out Martins 2005b AA within tolerance in all five runs (RMS 0.07-0.27; A 0.13-0.29).
- Acetic-acid volatility is not a cause: Kaw ~1e-5 at 25 C, ~1e-3 at 120 C (textbook Sander value, not in
  henry_constants.yml); < 0.5 % in the gas of a sealed vial.
- But arm A still identifies log10 k_aa_frag -2.60 +- 0.06 at lower total cost, driven by non-AA rows (Davidek
  2002 FA + AA, Martins 2003b FA).

**Verdict: EXPLAINED for the acetic-acid data (M1: the loss stands in for pH-driven slowing in a fixed-pH trunk),
P 0.6; the constant's residual identification comes from another defect, so the register stays OPEN until arm E.**
Consequence: when the trunk gets a pH state, drop or re-test k_aa_frag first.

## 5. CF-012 residual

The source conflict stays explained. The residual is not a pH-drift artefact: at in-situ pH the engine's Hofmann
contrast is 3.71 vs 5.27 -> 4.47 (not 3 vs 5), so its 55x norfuraneol pH-5 / pH-3 ratio implies an even steeper NF
pH law (~1.4 dex per in-situ unit) than the label reading, against Meynier's NF falling with pH. The engine's pH-5
Hofmann pot probably also runs too acidic (it fails Meynier's negative control by the same CBX mechanism). Residual
stays MOD (NF pH law + amine fate), not pH drift. Hou's trajectory (fast collapse in 20 min) says the Zhou / Jiangnan
pots spend almost the whole hold below pH 5, which the engine misses by ~1.3 units at 20 min.

## 6. Unbuffered pots without a measured pH record (ledger scan, heuristic on `medium` / `ph`)

Unbuffered, >= 80 C: 297 entries; 254 have no pH record (12 in_model_fit, 14 in_model_check, 180 digitised, 28
figure_only, 16 text_only, 4 supplement_missing); 43 have a final or running pH. Aqueous model pots in fits or checks
with no pH record: zhai2019_fig3a_ttca_cama (fit), zhai2019_fig3b (check), zhai2021_fig3a_ttca_tladder (fit),
zhai2023foodchem_table1_ttca_tladder (fit; only an unattached "7.0 -> 4.9" in text), kang2026 S3 / S4 / S4b (fit;
the endpoint is in kang2026_text_final_ph, not linked), zhou2023c_fig1 (check), deng2022_table1 (fit),
zhou2024_table2 (fit; a pH figure exists, zhou2024_fig3_fig4, figure_only: digitise), zhang2020_table1 (fit),
xia2022_fig5 (check), yu2018_table1 (check, pH 10), lin2021 / lin2022 text rows (check), schirle1992 (fit; pH
measured, not printed). Not yet in the ledger: wei2022 / zhou2022 / zhou2024b / cui2022 GX- and AX-ARP pots (same
design as Hou 2026, which measured 7 -> 5.05 in 20 min) and Hou 2026 itself.

## 7. A minimal declared drift model

1. Per pot: label pH and its reading temperature, buffer kind and molarity, titratable charge at MEASURED constants
   (cysteine GKL, ARP at the DFG analogue 8.18 / 2.20 declared, TTCA parent acid, NH4+ and H2S when charged) -> SID
   (exists).
2. One latent per time: net strong-acid equivalents n_H(t). Where a cooled pH is measured, invert the cooled charge
   balance for n_H(t_i) (no kinetics), interpolate monotonically, and solve the in-situ pH at T from the same
   inventory. Report the formic / acetic composition band (+-0.4 pH).
3. Where no record exists: n_H(t) from a class law (acid equivalents per mol precursor consumed, with a lag), its
   band set by leave-one-pot-out over same-class pots with records: pentose-ARP in water at 120 C (Hou, Zhou 2023,
   Zhang 2025, Zhou 2024 once digitised, Kang final); sugar + amino acid in phosphate (Martins 2005b, Cao, Deng 2021,
   Nguyen); sugar + cysteine in pyrophosphate / phosphate (Meynier <= 0.2 as a bound).
4. A vessel operator for H2S (and NH3 if wanted): dissolved fraction 1 / (1 + alpha0 Kaw(T) Vg/Vl).
5. Pass rules before use: cooled-pH RMS <= 0.3 leave-one-pot-out within each class; Meynier within 0.2; the
   operator, not the engine's acid kinetics, feeds the rate laws until the acid / amine budget passes both.
