# Pyrazine rivals: verdicts (2026-10-31)

Files in this directory:
- `prereg.yaml`: written 13:48-14:20, plus amendment_1 at 14:30. Both came before any run. The prereg is NOT blind, and its disclosure section says which values were already known.
- `pyrazine_rivals.py`: the power checks and the scoring. Seed 20261031, 4000 draws (20000 for S3).
- `pyrazine_rivals_result.{json,md}`: the results.
- `posthoc_summary.{py,json}`: the per-system minimum z (the lab-level rule in the prereg), and one POST HOC diagnostic for P1x.
- `run1.log`: the run log.

No file under `src/` or `data/` was edited, and nothing was committed.

## Power and size, then scores

| test | lab(s) | power (material) / size | powered? | result |
|---|---|---|---|---|
| S1: Leahy mixing bound, R = MPZ²/(PZ·DMP) ≤ R_max | Reineccius | 0.86-0.93 against P1x; 0.82-0.89 against P2; size 0.02 | yes | 14 of 15 independent rows are above the bound. Lys-Glc R = 11.7 / 14.4 / 17.3 against R_max 4.16, z = 3.6 / 4.1 / 4.5. Per-system minimum z: Lys-Glc 3.6, Lys-Rib 7.4, Asn-Glc 3.2, Asn-Fru 2.2, Lys-Fru 0.65. At the top of the RF band: 2.3 / 4.9 / 2.0 / 1.5 / 0.4. |
| S2: same channel, second lab | Ho network | ARP pot rows have power 0.54-0.70 (DMP25 is below the calibration intercept), so they are dropped. The Cys, Cys+MGO and GSH rows are powered. | partly | No powered Ho row rejects (maximum z 1.7, for Cys+MGO at pH 6 and 90 min). Every Cys row has R < 4, which no rival forbids. |
| S2u: uncalibrated sets | Hwang, Martin, Deng, Van Lancker, Shu | power ≤ 0.50 under the ±0.7 dex area-bias band | no | Reported only. Hwang's Lys pot (R = 94) is still above 4 at the worst bias (z 2.7). |
| S3: 15N identity, D = P2/f_A + P0/(1-f_A) - 1 | Ho (Zhou 2024) | 0.97 (primary band); 0.995 (secondary); size 0.10 | yes | D_obs = 0.129 / 0.093 / 0.076 / 0.076. **Primary band ([M-H]⁺ leakage ≤ 0.30): not rejected.** D_crit is 0.139; leakage at the band top alone gives D = 0.149. **Secondary band (correction assumed): rejected at all four times** (D_crit 0.057). |
| S4: carbon origin | Nestlé, Reading | NC 1.00; P2a_aa 0.99 | yes | MPZ and DMP25 are 0 % labelled (Amrani) and ≤ 7 % (Low). NC fails, as a negative control should. P2a_aa is limited to ≤ 15-30 % of MPZ, while S1 needs ≥ 42 %. |
| S5a: pyrrole mass balance | Ghent | 0.95 / size 0 | yes | Pyrrole / (2·DMP) is 19-116 in all five pH 4 rows (B = 1.28-2.06 dex, about 6-9 σ). The pH 7 rows do not reject. |
| S5b: sign of pyrrole/DMP against MGO:Ala | Ghent | (band-free) | yes | 38 → 55 → 231, rising, as P1, P1x and P2 predict. |
| S6: Strecker yield vs the P0 constant | DFA vs Ho | 0.38 | **no, dropped** | MGO: P0 under-predicts by 1.46 dex at the most favourable corner of its bands, and by 2.33 dex at defaults. GO: 0.29 / 1.19 dex. |
| S7: pH 4 vs 7 DMP, fed MGO | Ghent | 0.93 / 0.04 | yes | Observed +0.32 dex against P0's −1.74. P0's pH placement is rejected (2.06 dex against a 2σ of 0.73; the acetate vs phosphate buffer difference is inside that band). |

**Post hoc, not scored.** Under P1x, one value of u would have to serve all sugars at a given T, because the aminoketones are the same. Lys-Rib then forces u ≥ 46-68. That in turn needs h = 0.09-0.36 in the Lys-Glc, Lys-Fru, Asn-Glc and Asn-Fru pots, which is outside the pre-registered h_linear band [0.5, 1]. This needs its own prereg, but it points against P1x.

## Verdicts (rule: admitted if every powered test passes within the bands, provisional if single-lab; rejected if a powered test fails in ≥ 2 labs, or in 1 lab plus a mass balance)

- **P0 (engine): REJECTED.**
  - Ghent mass balance (S5a) and Ghent's within-study pH contrast (S7).
  - Reineccius S1, because the engine uses one shared k_cond.
  - S6 points the same way but is underpowered.
- **P1 (one shared condensation constant): INCONCLUSIVE by rule.**
  - It fails S1 in one lab, and S1 is a constant-free Cauchy-Schwarz bound, so only a response-factor bias of 0.4-1.1 dex could rescue it.
  - No powered second-lab row exists.
  - P(P1 is wrong) ≈ 0.85.
- **P1x (pair-specific constants): ADMITTED, PROVISIONAL (single lab), fragile.**
  - It is rejected if Zhou 2024 corrected its isotopomers for [M-H]⁺ (S3 secondary).
  - The post hoc cross-sugar diagnostic also argues against it.
- **P2b (a sugar-carbon MPZ route that does not go through A × M pairing; for example a hydroxy-C3 aminoketone with dehydration, or a three-component route): ADMITTED, PROVISIONAL.** It passes every powered test. Its discriminating channels (S1, S3) are single-lab each.
- **P2a_aa (the amino acid's Strecker aldehyde adding to the A+A dihydropyrazine to give MPZ): REJECTED** in the 180 C matrices: two labs, same channel (S4). Water was not tested.
- **P2a_sugar: inconclusive.** It is not separable by label. It is disfavoured because aldehyde addition to the bare dihydropyrazine stays ≤ 15-30 % of MPZ in pots where HCHO addition demonstrably happens.
- **NC: fails S4** (control works).

## Missing experiments (cheapest first)

1. **Zhou 2024's isotopomer method.** Was [M-H]⁺ corrected? Read the SI (10.1021/acs.jafc.4c03706.s001, a Crossref component) or ask the authors. This alone decides between P1x and P2. It is a reading task, not bench work.
2. **Fed GO + MGO together** at 2-3 ratios with Ala, plus one 15N-Ala arm and one [13C3]MGO arm. Measure PZ, MPZ and DMP25/26 with SIDA or calibration, plus aminoacetone and aminoacetaldehyde by derivatised LC-MS/MS.
   - P1x predicts R independent of the ratio.
   - P2b needs a sugar-derived C3 and predicts R ≈ 4 in a sugar-free pot.
   - This extends EX design E5 in EVIDENCE_MAP sec. 3.
3. **Fed hydroxyacetone / DHA + GO + amine, argon vs air.** This tests the dehydration route (P2b): MPZ rises under argon while PZ and DMP fall.
4. **A cross-lab replicate of Leahy's Lys-Glc pH 9 pot** with in-matrix calibration (this removes the single-lab status of S1).

## Papers (DOIs Crossref-verified this session)

- Zhou et al. 2024, JAFC 72:18630, 10.1021/acs.jafc.4c03706; SI 10.1021/acs.jafc.4c03706.s001 (the isotopomer method).
- Zhou et al. 2024a, JAFC 72:5878-5886, 10.1021/acs.jafc.4c00021 (methods ref. 20 of Zhou 2024).
- Shibamoto & Bernhard 1977, Agric. Biol. Chem. 41:143-153, 10.1271/bbb1961.41.143 (dehydration ring closure).
- Koehler & Odell 1970, JAFC 18:895-898, 10.1021/jf60171a041.
- Low et al. 2006, JAFC 54:5976-5983, 10.1021/jf060328x (absolute concentrations for low2007).

Lead (not re-verified here): Hwang, Hartman, Rosen & Ho 1993, the ion-cluster isotopomer method, whose DOI was not retrieved (rate-limited). For the other pyrazine leads, the review lists Shibamoto & Bernhard 1976, Scalone 2015, Van Lancker 2010 and Hayashi & Namiki 1980 as Crossref-checked by an earlier pass.
