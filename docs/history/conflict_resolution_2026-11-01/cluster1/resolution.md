# Cluster 1: same-lab reprints and text / table slips (CF-018, 019, 020, 029, 040, 057, 067, 070, 074)

Done 2026-11-01. Read-only on tracked files. Evidence files: `scratchpad/c1/` (pdftotext layers of hofmann1999 / 2000 /
2000b, rainercremer2000, weenen2001, Xia2022, devleeschouwer2008 / 2009b and Zhai2023; `fingerprint_scan.py`,
`value_only_scan.py`, `fingerprint_hits.json`, `fingerprint_hits2.json`, `value_only_hits.txt`). Also
`scratchpad/cf018/` and `scratchpad/cf029/`, the PDF text layers from the two sub-agents that checked CF-018 and CF-029
cell by cell. I spot-checked their key claims; for example, the "pairs" in CF-018 that are one PDF do name the same PDF
and DOI. Quotes are under 15 words. P is my probability that the stated cause is right.

**The common decision rule, stated before reading the evidence:**
- **DUPLICATE**: at least 3 numeric cells (value, and SD where printed) agree to the printed precision. A single
  stated unit factor common to every cell is allowed.
- **NOT A DUPLICATE**: at least 3 shared-condition cells differ beyond rounding or read error.
- **UNDECIDED**: fewer than 3 comparable cells exist.

For self-contradictions:
- **SLIP**: the two numbers differ by a unit factor, a power-of-ten / decimal shift or a digit transposition, and the
  paper's own arithmetic (totals, derived columns, a second print) picks one.
- **RUN DIFFERENCE**: the discrepancy cannot be produced by any single slip, for example opposite directions in two
  cells, or a constant non-unit factor.

---

## Summary

| CF | verdict | P | proposed status |
|---|---|---|---|
| 018 | 4 of 7 "pairs" are two dossiers of ONE PDF (repo-side duplicate dossiers); edwards / wedzicha ME + GSH rows and the cammerer 90 C row are genuine reprints; martins2003 / 2003b undecided (figure-only) | 0.97-0.99; martins 0.97 untested | duplicate_confirmed (martins pair stays pending) |
| 019 | Leahy ch. 7 / ch. 18: one data set (12 k values + n identical). zhou2024 / 2025b, mikami2017 / nakamura2020, lee2022 / lee2024: not duplicates | 0.99 / 0.95 / 0.9 / 0.9 | duplicate_confirmed for Leahy; the rest not_an_error (split) |
| 020 | xia2022 RESOLVED (zero-order slopes; "pseudo-first-order" is the slip); cremer2000 EXPLAINED (oven range, not the isotherm set); ge1997 R2 item RESOLVED (threshold is on r); weenen2001 run difference, not a slip; kim1996 OPEN; agcam2022 explained by its dossier; hofmann2000 already ruled | see sec. 3 | split by paper |
| 029 | Selke / Frankel autox column, Crowther T7 / Aspelund, Snel / Snel thesis: duplicates. Hidalgo / Zamora, Xin 2026b / Xin 2026: NOT duplicates. B6 does not read Selke; only one cell is read once (lipid.py:689) | 0.97-0.99; 0.95-0.97 | duplicate_confirmed (3 pairs); strike 2 pairs |
| 040 | DV 2009 Part I glucose column = DV 2008 a_w 0.92 column, every cell; the "different k units" is a typo in Part I ("mm-1" for "min-1"). DV 2009 Part II's fixed (c) cells are a third print | 0.99 | duplicate_confirmed |
| 057 | RESOLVED: the value is 26.000 (trailing zeros dropped); the class subtotal reproduces it exactly | 0.97 | not_an_error |
| 067 | Figure / table: the CD-factor / transcription-slip explanation is REJECTED (constant 0.88 factor, not a power of 2); a different analysis set survives. OPD: OPEN; the "2 mL; 1 mmol/L" text recurs in hofmann2000, so it is the lab's written protocol (misprint P 0.5) | 0.6 / 0.5 | open (explanations narrowed) |
| 070 | One phenylalanine series printed twice (3 of 5 cells identical to the printed digit; equal loss from MGO and glucose is chemically implausible). Probably measured in the glucose pot | 0.9 duplicate; 0.7 owner = 2000b | duplicate_confirmed |
| 074 | Only 2 cells, so UNDECIDED by the cell rule; likelihood gives P(one run) ~ 0.8. Control: ARP-Phe air / argon PA ratio 9.2 is 0.42 dex from 3.5, beyond 0.15, so B29's 3.5 is flagged as a possible vessel comparison | 0.8 | explained (flag kept) |

---

## 1. CF-018: trunk reprints (checked cell by cell by a sub-agent; spot-checked here)

**Rule:** the common rule.

**Test:** the sub-agent lined up the dossier tables, ran pdftotext on the PDFs where needed, and grepped `src/` and
`scripts/` for readers. I re-checked the PDF / DOI identity of pairs 1-4 (grep of each dossier's "Source on disk" and
DOI).

| pair | test result | verdict |
|---|---|---|
| kocadagli2016nacl / kocadagli2016jafc | Both dossiers read `Kocada2016.pdf`, doi 10.1021/acs.jafc.6b01862. Not a reprint: one paper with two dossiers. The nacl dossier transcribes no table. | RESOLVED, repo-side duplicate dossier, P 0.99 |
| goncuoglu2016 / goncuoglutas2016 / goncuoglutas2017 | All three read `Goncouglu2016.pdf`, doi 10.1016/j.foodchem.2016.11.159. Table 1 cells agree in all three (k1 6.9 +- 0.8 / 15 +- 3.1 / 22 +- 2.0; k17 371 +- 241; ...). The fingerprint scan finds 27 shared value +- SD cells. | RESOLVED, one paper with three dossiers, P 0.99 |
| gursul2020 / gursulaktag2020 | Same PDF (`Gursul2020.pdf`, doi 10.1016/j.foodchem.2020.126620). Identical cells (k1 apple 37 C 123.7 +- 10.03); the scan finds 81 shared numbers. | RESOLVED, one paper with two dossiers, P 0.99 |
| huang2016 / huang2017 | Same PDF (`huang2016.pdf`, doi 10.1111/jfpe.12375). Table 1 identical (k1 5.68E-03 +- 6.74E-04 ...; 18 shared value +- SD cells). | RESOLVED, one paper with two dossiers, P 0.99 |
| edwards1992 Table 1 / wedzicha1991 Table 1 | ME and GSH rows identical in both PDFs (k1 0.029 / 0.020; k2 1.31 / 1.35; no SDs printed). The NAC and sulfite rows are new in 1992. | RESOLVED, partial reprint (2 rows), P 0.97 |
| cammerer1994 90 C row / wedzicha1992 | Cammerer footnote c credits "Wedzicha & Kaputo (1992)". a 1.06 = 1.06. Cammerer's "y/a" 2.96 is Wedzicha's y (Wedzicha's y/a is 2.79). The row maps to Wedzicha Table 5, not Table 3. | RESOLVED by the paper's attribution (only 2 cells), P 0.99 |
| martins2003 Figs 4-7 / martins2003b | Figure-only on both sides. Coarse reads agree within read error. Part II's methods hold only "2.1. Kinetic modelling". | OPEN, P(duplicate) 0.97. Decided by point-by-point digitisation, >= 3 matching (t, value) |

**Double-count check:**
- Pairs 1-4 are never read twice by code. Kocadagli: only the jafc entries (parameters_furanic / dicarbonyl). Hazelnut:
  anchored once on the goncuoglu2016 dossier (DISPUTED_SINK_BANDS). Juice: B7 hold-out H2, once. Huang: no reader.
- The double count is in the ledger tally only: three in_model_check hazelnut entries and two juice entries.

**Cause:** the four "same-lab reprints" are a repo-side TU error (one PDF dossiered under two or three stems), not a
paper-side reprint. This should be recorded so the same check is run at dossier-creation time.

```yaml
  - id: CF-018
    status: duplicate_confirmed      # martins2003 pair excepted (pending digitisation)
    sources: [kocadagli2016jafc, kocadagli2016nacl, goncuoglutas2017, goncuoglu2016, goncuoglutas2016, gursulaktag2020, gursul2020, huang2017, huang2016, wedzicha1991, edwards1992, wedzicha1992, cammerer1994, martins2003b, martins2003]
    explanations:
      - "Pairs kocadagli, goncuoglu (x3), gursul, huang: one PDF dossiered under two or three stems (repo-side), not a re-publication: same file and DOI."
      - "edwards1992 ME / GSH rows and cammerer1994's 90 C row: paper-side reprints of wedzicha1991 Table 1 and wedzicha1992 Table 5."
    tests:
      - "2026-11-01 cell check against the PDFs: identical DOI / file for the four one-PDF pairs; ME / GSH k1, k2 identical (4 cells); Cammerer footnote c credits Wedzicha & Kaputo 1992; a 1.06 = 1.06, y 2.96 = 2.96."
      - "Code grep: no fit or check reads two copies of any pair (kocadagli jafc only; hazelnut anchored once on goncuoglu2016; juice once in B7 H2)."
    remedy: >-
      copy_of (fitting.count_once): kocadagli2016nacl_fig3_fig4_time_courses -> kocadagli2016jafc fig3 / fig4 entry (the
      nacl dossier has no Table 1 / 2 copy; drop kocadagli2016jafc_tab1_tab2 from this pairing);
      goncuoglu2016_hazelnut_roasting_matrix, goncuoglutas2016_hazelnut_roasting_matrix -> goncuoglutas2017_hazelnut_roasting_matrix;
      gursul2020_fruit_juice_storage_matrix -> gursulaktag2020_fruit_juice_storage_matrix;
      huang2016_tab1_tab2_leu_ile_methylbutanal_constants -> huang2017_tab1_tab2_leu_ile_methylbutanal_constants;
      edwards1992_table1_3dg_activation_four_nucleophiles (ME, GSH rows only) -> wedzicha1991_table1_3dg_activation_glycine_catalysis;
      cammerer1994_table3_melanoidin_composition_conditions (90 C row) -> wedzicha1992 Table 5 (1.0 : 1.0 row; Cammerer
      prints y in its y/a column and pH 5 for 5.5). Retire the duplicate dossiers (kocadagli2016nacl, goncuoglu2016 or
      goncuoglutas2016, gursul2020, huang2016) or mark them `same_pdf_as:`; one ledger entry per data set.
      martins2003_figs4_7 -> martins2003b: pending point-by-point digitisation (no fit reads Figs. 4-7).
```

## 2. CF-019: trunk same-lab pairs

**Rule:** the common rule. For CF-019 the "original" is the earlier publication.

### 2a. leahy1989 (ch. 7, Feb 1989) / leahy1989a (ch. 18, Oct 1989): RESOLVED duplicate, P 0.99

**Test:** cell comparison of the lysine-glucose pH 9.0 block in the two dossiers' Table I. ch. 18's Table I was
verified cell by cell against a raster.

- All 12 k values agree to the printed digit, and so do all n. Pyrazine 95 / 85 / 75 C: 3.596 / 0.490 / 0.214;
  methylpyrazine: 2.837 / 0.422 / 0.159; 2,5-DMP: 0.186 / 0.0247 / 0.00668; 2,3-DMP: 0.0229 / 0.00309 / 0.000677.
- Five secondary cells differ: an r2 (0.995 vs 0.976), three 2,3-DMP intercepts (-0.00173 vs -0.00017, which looks
  like dropped digits), and two Ea (0.1 and 0.4 kcal/mol). A refit of the identical k reproduces ch. 7's Ea.
- These are typescript slips in a reprint, not a second run.
- The ch. 7 Table III 2 h total, 13.1 ppm, is the same in both.

**Double-count check:**
- The B18 fit reads ch. 18's within-study ratios k(pH x) / k(pH 9) (`generate_kinetic_core_b18_fit.py` LEAHY_RATIOS).
- The B18 ship rule T3 / T4 and B21 T5 read ch. 7 Table III (the pH 9 distribution and total).
- No objective reads both copies. The pH 9 arm is the denominator of a fitted ratio AND the subject of a ship check
  (a different observable of the same experiment). Record that the check is not independent of the fit.
- Ledger: both Table I entries are `in_model_fit`. Nothing fits ch. 7 Table I.

### 2b. zhou2024 / zhou2025b: RESOLVED not a duplicate, P 0.95

**Test:** CF-060's adjudicated cells, read from page vectors and a 600 dpi raster.
- Every dicarbonyl differs 7-24x at matched times: 3-DX 3.7-5.3 vs 0.37-0.52 mM; GO 1.0 -> 4.9 vs 0.15 -> 0.31.
- 2,6-DMP is detected in 2025b and "not detected" in 2024 (zhou2025b dossier flag 7).
- The two papers share the design only (the ledger note already says so).

### 2c. mikami2017 Fig. 6B / nakamura2020 Fig. 3b: RESOLVED not a duplicate (an independent same-lab replicate), P 0.9

**Test:** both dossiers digitised the NF-alone pot at pH 8, 100 C, with OPD added after heating.

| t (h) | NF Mikami | NF Nakamura | DA Mikami | DA Nakamura | MGO Mikami | MGO Nakamura |
|---|---|---|---|---|---|---|
| 0.5 | 0.63 | 0.735 | 0.030 | 0.028 | 0.011 | 0.038 |
| 1 | 0.25 | 0.205 | 0.063 | 0.052 | 0.004 | 0.018 |
| 6 | ~0 | ~0 | 0.046 | 0.061 | ~0.001 | 0.002 |

- At least 5 cells differ beyond read error (+-0.02-0.03 mM NF; +-0.003 mM quinoxalines).
- The trapping step differs: OPD "5 h at 20 C" (Mikami) vs "1 h RT" (Nakamura).
- The decay constants agree: 0.032 min-1 (Mikami) vs 2.062 h-1 = 0.034 min-1 (Nakamura 3b), within 7 %.

So this is a same-lab replicate that supports the B61 anchor; it is not a copy. Counting it as a second
INDEPENDENT-LAB point would be wrong (one lab, one protocol).

### 2d. lee2022 / lee2024: RESOLVED not a duplicate, P 0.9

The same model-cake experiment published in two halves with disjoint species: non-volatile markers vs oven-vapour
volatiles. The lee2024 dossier: "the same matrix, oven and formulations ... one dataset". There are no cells to
double-count. Shared conditions only, so carry `same_experiment_as`.

```yaml
  - id: CF-019
    status: duplicate_confirmed     # the Leahy pair; the other three pairs: not_an_error, moved to a note
    sources: [leahy1989, leahy1989a]
    explanations:
      - "leahy1989a's pH 9.0 rows (Table I) and pH 9.0 column (Table III) reprint ch. 7's lysine-glucose block."
    tests:
      - "2026-11-01: 12 / 12 k values and n identical to the printed digit; 5 secondary cells (r2, 3 intercepts, 2 Ea) differ by dropped digits / 0.1-0.4 kcal/mol (refit of the identical k gives ch. 7's Ea)."
      - "Not duplicates (same cells differ beyond read error): zhou2025b vs zhou2024 dicarbonyls 7-24x (CF-060); nakamura2020 Fig. 3b vs mikami2017 Fig. 6B NF 0.735 vs 0.63 mM, MGO 0.038 vs 0.011 mM at 0.5 h, different OPD step; lee2024 vs lee2022 disjoint species."
    remedy: >-
      copy_of: leahy1989a_tab1_lys_glc_pyrazine_rates_ph_ladder (pH 9.0 rows) -> leahy1989_tab1_lys_glc_pyrazine_rates
      (lysine-glucose rows); leahy1989a Table III pH 9.0 column -> leahy1989 Table III. B18 reads only ch. 18's ratios;
      B18 T3 / T4 and B21 T5 read ch. 7 Table III: same experiment, different observable -- record that these ship
      checks are not independent of the fit. Change leahy1989_tab1 status to in_model_check (nothing fits ch. 7 Table I).
      Not duplicates: zhou2025b (same design, not a re-plot; CF-060); nakamura2020 Fig. 3b (a same-lab replicate of
      mikami2017 Fig. 6B: k 0.034 vs 0.032 min-1; one lab, never two); lee2024 (same experiment, disjoint species).
```

## 3. CF-020: self-contradictions in trunk papers, split by paper

| paper | rule outcome and test | verdict |
|---|---|---|
| hofmann2000 | Already ruled (0.5 mmol/L is a misprint; it cannot hold pH 7). | RESOLVED earlier |
| **xia2022** | Rule: the unit is right if the axis can hold the rate. Fig. 5B (3-DG) is drawn to 15 mmol/L over 0-80 min. As a first-order constant on 200 mM glucose, k(3-DG, 130 C) = 0.0357 min-1 would give ~190 mmol/L by 80 min. As a zero-order slope, 0.0357 mmol L-1 min-1 x 80 min = 2.9 mmol/L, inside the axis. G-ARP at 130 C: 0.43 x 35 min = 15 mmol/L, plausible. As min-1 it would mean near-total conversion of 200 mM in minutes. | **RESOLVED**, P 0.9. Table 1 / 2 k are zero-order initial slopes (mmol L-1 min-1, at fixed 200 + 200 mM); "pseudo-first-order" in the prose is the slip. |
| **cremer2000** | The "70 to 110 C" sentence reads "heated (thermostatted) for 3 to 120 min in the headspace oven at 70 to 110 C". It is the oven's range for ALL experiments, matched by the autosampler line "sample temperature, 70-110 C; thermostatting time, 3-120 min". Plant powders ran at 90 C (Fig. 3 caption); the isotherms are 80 / 90 / 100 / 110 (Fig. 6), with four Arrhenius points (Fig. 7). | **EXPLAINED**, P 0.85. Not a contradiction about the isotherms; there is no 70 C series. |
| **ge1997** (R2 item) | The stated threshold "rs = 0.714" is a critical value of a correlation coefficient. Table 2 prints R2. sqrt(0.6046) = 0.778 > 0.714; every printed R2 (0.6046-0.7417) gives r 0.78-0.86 > 0.714. | **RESOLVED** for this item, P 0.8: the text compares r, the table prints R2. The K vs k1 inconsistency at 68 / 81 C (implied 1.19 / 1.13 vs 0.71 / 0.96: not unit factors or transpositions) stays **OPEN**. |
| **weenen2001** | Rule: a slip if one factor maps one table onto the other. Both tables nominally use one protocol: Likens-Nickerson, citrate 50 mM, ~1.3 mmol carbonyl + Phe, 2 h. Glucose + Phe is 0.006 % (Table II) vs 0.03 % (Table III), x5. Glc-Phe ARP alone is 1.6 % (Table II) vs 0.5 % (Table III), x0.31. Opposite directions, so no single slip maps one table to the other. The text supports Table II's glucose cell ("degradation of a phenylalanine solution without glucose results in the same amount"; Phe alone 0.006 %). The Methods say pH 3.5, the tables pH 3.2. | **RUN DIFFERENCE, OPEN**, P(single slip) 0.15. Any use of the ARP / glucose contrast must carry both: 270x (Table II) vs 17x (Table III), a 1.2 dex band. The evidence map K2 currently quotes 270x only. |
| kim1996 | 1 % (Fig. 4 caption) vs 10 % (text) glucose. No internal arithmetic on disk decides it. | **OPEN**, P 0.5 either way. Decided by Fig. 4's M-1 level vs the 1 % spikes of Figs. 2-3 (same matrix): a 10x loading should show. |
| agcam2022 | Methods 0-90 min vs prose 0-120 min. Three t1/2 inconsistent with their own k under every unit reading (dossier sec. 3). | **EXPLAINED** by the dossier's rulings (methods authoritative; t1/2 refused), P 0.85. |

```yaml
  - id: CF-020
    status: open      # proposed split: xia2022 / ge1997-R2 / cremer2000 / agcam2022 closed, weenen2001 / kim1996 / ge1997-K open
    tests:
      - "xia2022 (2026-11-01): Fig. 5B axis 0-15 mmol/L over 0-80 min; k(3-DG, 130 C) 0.0357 as min-1 on 200 mM would reach ~190 mmol/L, as mmol/L/min 2.9 mmol/L: the printed unit is right, 'pseudo-first-order' is the slip."
      - "cremer2000: '70 to 110 C' is the headspace oven's range for all experiments (paired with '3 to 120 min'; powders at 90 C); isotherms 80-110 C (Figs. 6-7)."
      - "ge1997: the 0.714 threshold is on r; sqrt(0.6046) = 0.778 > 0.714 (every R2 in Table 2 passes as r)."
      - "weenen2001: Table II vs Table III differ x5 (glucose) and x0.31 (ARP) in opposite directions: two runs, not one slip; the text supports Table II's glucose = Phe-alone 0.006 %."
    remedy: >-
      xia2022: rows in mmol L-1 min-1 as zero-order initial slopes at fixed 200 + 200 mmol/L. cremer2000: isotherms 80 / 90 /
      100 / 110 C. ge1997: R2 item closed; K / k1 at 68 and 81 C remain condition_uncertain. weenen2001: both tables kept as
      separate runs; the ARP / glucose PA ratio at pH 3.2 carried as the band 17-270x (evidence map K2 to be corrected).
      kim1996, agcam2022: as before.
```

## 4. CF-029: lipid-lane reprints (checked against the PDFs by a sub-agent)

**Rule:** the common rule. The verdict is taken per column where a table mixes copied and new columns.

| pair | test result | verdict |
|---|---|---|
| selke1978 Table I / frankel1981 Table II, autoxidation column | All 20 cells identical (heptane 4.4, octane 2.7, octanal 11, nonanal 15, ...). The footnote: "Data from ref. 21." The photosensitized column is independent. | RESOLVED duplicate (column), P 0.99 |
| aspelund1983 Table I / crowther1980 Table 7 "Aspelund, 1978" column | 6 / 7 -dH identical (6.04, 8.11, 8.89, 11.25, 13.89, 18.06). The 7th (pentanone 2.92) is a thesis value. | RESOLVED duplicate (column), P 0.98 |
| snel2023 / snel2023b (thesis ch. 6) | The thesis says "This chapter has been published as Snel ... LWT". 40 / 40 K values identical with SDs. The Ratio column is recomputed (10 cells off by 1-3 units). | RESOLVED duplicate, P 0.97 |
| zamora2015 Table 1 / hidalgo2016 control | PEA N2 30.99 +- 3.37 vs 21.9 +- 3.6; PEA air 11.23 vs 10.6; PAC air 126.82 vs 129.7. 3 / 4 means and all SDs differ. | RESOLVED **not** a duplicate (same-lab re-run), P 0.95 |
| xin2026 control / Xin2026b ENPI arm | Totals 2433.19 vs 193.95 (aldehydes), 6778.15 vs 3079.28 (ketones); 10-23x apart. This pair IS the evidence for HS_SPME_SAME_SAMPLE_DISPERSION (10, 23). | RESOLVED **not** a duplicate, P 0.97 |

**B6 / B28:**
- B6 fits only Frankel 1989 (FRANKEL_ZERO_ADDITIVE). Selke 1978 is read nowhere in `src/` or `scripts/`.
- FRANKEL1981_OLEATE_SLATE is read once (`lipid.py:689`, `_b28_yields`: the autoxidised NONANAL 15 as a fixed share).
- There is no double count in any objective, only two `in_model_fit` ledger rows for one measurement.
- The register text "the Selke slate in B6" is wrong.

The sub-agent also found a dossier slip: `xin2026_extraction.md` l. 59 gives 922.04 ug/kg as hexanal, but the PDF
gives it as benzaldehyde.

```yaml
  - id: CF-029
    status: duplicate_confirmed      # three pairs; two pairs struck
    sources: [selke1978, frankel1981, aspelund1983, crowther1980, snel2023, snel2023b]
    observation: >-
      (rewrite) Frankel 1981 Table II autoxidation column = Selke 1978 Table I (20 / 20 cells; footnote 'Data from ref.
      21'); crowther1980 Table 7 Aspelund column = aspelund1983 Table I (6 / 7 cells); snel2023b ch. 6 = snel2023 (40 / 40
      K with SD). B6 does not read Selke; B28 reads one Frankel cell once (lipid.py:689).
    tests:
      - "PDF cell checks 2026-11-01; grep: selke1978 read nowhere, FRANKEL1981_OLEATE_SLATE read once."
      - "Not duplicates: hidalgo2016 control vs zamora2015 (PEA N2 21.9 vs 30.99, PAC air 129.7 vs 126.82; all SDs differ); Xin2026b ENPI vs xin2026 control (10-23x apart; the HS-SPME dispersion evidence)."
    remedy: >-
      copy_of: frankel1981_tab2_oleate_lOOH_volatile_slates (autoxidation column only) -> selke1978_tab1_oleate_lOOH_slate;
      crowther1980_tab3_4_igc_processed_soy (Table 7 Aspelund column) -> aspelund1983_tab1_3_igc_soy_80_100c;
      snel2023b_fig6_2_raw_rhc_aldehydes -> snel2023_fig2_3_rhc_vs_protein (the thesis digitisation stays usable as the
      original's values). Remove selke1978's reverse duplicate_of flag (circular). hidalgo2016 / zamora2015: same_lab_rerun,
      no copy_of. Xin2026b / xin2026: keep same_material_as, no copy_of. Fix xin2026 dossier l. 59 (922.04 is benzaldehyde).
```

## 5. CF-040: De Vleeschouwer 2008 / 2009 Part I (and Part II)

**Rule:** the common rule, plus a unit check. If the units differ, every cell must differ by the same factor; if the
cells are identical, the "unit difference" is a label slip.

**Test:** pdftotext of both PDFs, Table 2 (2008) vs Table 3 (2009b), glucose column.

| cell | DV 2008 a_w 0.92 | DV 2009b glucose |
|---|---|---|
| k_Fref | 3.57 +- 1.38 | 3.57 +- 1.38 |
| k_Eref | 0.10 +- 0.04 | 0.10 +- 0.04 |
| k_INTref | 1.70 +- 1.05 | 1.70 +- 1.05 |
| k_Mref | 1.23 +- 0.49 | 1.23 +- 0.49 |
| k_Bref | 3.90 +- 3.68 | 3.90 +- 3.68 |
| k_Aspref | 26.43 +- 5.76 | 26.43 +- 5.76 |
| Ea_F | 159.2 +- 29.5 | 159.2 +- 29.5 |
| Ea_E | 113.2 +- 32.3 | 113.2 +- 32.3 |
| Ea_INT | 117.5 +- 25.2 | 117.5 +- 25.2 |
| Ea_M | 105.7 +- 29.1 | 105.7 +- 29.1 |
| Ea_B | 180.3 +- 38.5 | 180.3 +- 38.5 |
| Ea_C | -6.7 +- 0.2 | -6.7 +- 0.2 |
| Ea_Asp | 105.4 +- 10.6 | 105.4 +- 10.6 |
| Ea_X | 668.9 +- 35.2 | 668.9 +- 35.2 |

The fingerprint scan finds 12 shared value +- SD pairs independently.

**Units:**
- 2008 prints "kFref (10-3 min-1)".
- The 2009b text layer prints "kFref (10^-3 mm^-1)", while every other row of the same table reads "min^-1".
- Eq. 7 in both papers is first order in Int1 (d[AA]/dt = kF[Int1] - kE[AA]), so the unit is min-1.
- The "mM-1 min-1" in the 2009b dossier is a misreading of the typesetting typo "mm-1". The registry already carries
  `printed_unit_is_a_typo_10e-3_mm-1_means_min-1`, and the Part II dossier confirms it.

**Third print:** DV 2009 Part II Table 3 carries the same control cells (3.57, 0.10, 1.70, 1.23, 3.90, 26.43, 117.5,
...), each footnoted c ("fixed at their estimated values"). Part II is in_model_check.

**Double count:**
- B3 reads 2009b (k_E, Ea_E rows; k_F, k_INT, k_Asp operative via `_DV1_SOURCE`).
- B14 reads 2008's four columns only as ratios to the 0.92 column (the a_w band). The shared column is the
  denominator there, identically 1, so it carries no independent weight.
- No objective counts it twice. The ledger shows two `in_model_fit` entries.

**Verdict:** RESOLVED duplicate, P 0.99.

```yaml
  - id: CF-040
    status: duplicate_confirmed
    sources: [devleeschouwer2008, devleeschouwer2009b, devleeschouwer2009]
    tests:
      - "2026-11-01, both PDFs: all 14 value +- HPD cells of the glucose column equal DV 2008 Table 2 a_w 0.92 to the digit; 2009b's 'kFref (10^-3 mm^-1)' is a typesetting typo for min^-1 (every other row min^-1; Eq. 7 first order in Int1), so the units do not differ."
      - "DV 2009 Part II Table 3: the same cells, footnote c (fixed), a third print."
    remedy: >-
      copy_of: devleeschouwer2009b_tab3_glucose_acrylamide_rows -> devleeschouwer2008_tab2_aw_multiresponse_constants (a_w
      0.92 column); Part II's footnote-c cells -> the same. Drop the ledger flag unit_conflict_between_dossiers and correct
      the 2009b dossier's k_F unit (10^-3 min^-1). No code change: B3 reads one copy, B14 uses the column only as its ratio
      denominator.
```

## 6. CF-057: Zhai 2023 Food Chem Table 1, furfural 140 C TTCA + Xyl "26 +- 0.725"

**Rule:** a dropped decimal if the paper's own subtotal requires 2.6 or 0.26; a real value if the subtotal requires 26.

**Test (internal arithmetic):** the oxygen-heterocycle subtotal for 140 C TTCA + Xyl prints 32.226 +- 0.812. The other
members of that class in the same column are:
- 2-methylfuran 2.578
- furan 2.319
- 1-(2-furanyl)-ethanone 0.867
- 2(5H)-furanone 0.462

32.226 - 6.226 = **26.000** exactly. The same subtotal arithmetic reproduces exactly in the neighbouring columns:
- 120 C TTCA + Xyl: 1.182 + 0.298 + 8.497 + 0.178 + 0.189 = 10.344, as printed;
- 140 C TTCA: 12.757, as printed.

The grand total also closes: 163.432 + 10.296 + 32.226 = 205.954, as printed.

The SD is not implausible: 0.725 / 26 = 2.8 %, while the same row's SDs are 2.6 % (3.381 +- 0.089), 3.4 % and
2.7 % (11.039 +- 0.302).

**Verdict:** RESOLVED. 26.000 with the trailing zeros dropped; the letter 'a' (row maximum) is consistent. P 0.97.

```yaml
  - id: CF-057
    status: not_an_error
    explanations:
      - "A real value printed without its trailing zeros (26.000)."
    tests:
      - "2026-11-01: the oxygen-heterocycle subtotal 32.226 minus the column's other four members (2.578 + 2.319 + 0.867 + 0.462) = 26.000; the same subtotal arithmetic closes exactly in the 120 C TTCA + Xyl and 140 C TTCA columns, and 163.432 + 10.296 + 32.226 = 205.954 (Total). SD 2.8 %, as the row's other cells (2.6-3.4 %)."
    remedy: "Furfural 140 C TTCA + Xyl = 26.000 +- 0.725 ug/L, as printed."
```

## 7. CF-067: Hofmann 1999, glycolaldehyde figure vs table, and the OPD amount

### 7a. Figure vs table

**Rule:** a Table 3 slip if the figure / table ratio is a power of 2 (CD factors are 1 + 1 dilution steps, 2^n) or a
single-cell digit error. A different analysis set if the ratio is a constant non-unit factor confined to one species.

**Test:**
- The glycolaldehyde figure / table ratios at 10 / 60 / 120 min, from the keyed rows, are 0.884 / 0.902 / 0.856: a
  mean of 0.88, spread +-0.02, within the figure's read error of 3-7 %.
- MGO on the same lower axis agrees (0.99-1.03), so the axis calibration is not the cause.
- 0.88 is not a power of 2, so a CD-factor slip is ruled out.
- Three cells shifted by one factor cannot be a digit slip.
- It is not the glycolaldehyde response factor: the printed Rf is 0.96, and Rf squared is 0.92.
- Glycolaldehyde is quantified as its O-ethyloxime against [13C4]-3-hydroxy-2-butanone ("4 mmol ethoxamine
  hydrochloride was added"), NOT with OPD. The OPD question (7b) therefore cannot touch the glycolaldehyde rows.

**Verdict:** explanation 2 (a slip in Table 3's CD factor) is REJECTED. Explanation 1 (the figure comes from another
analysis or replicate set, or carries an unstated species-specific correction) survives, P 0.6. The status stays open,
with the explanations narrowed. The keyed remedy (table value, with 0.066 dex as the error) stands.

### 7b. The OPD amount

**Rule:**
- misprint if the same laboratory, for the same SIDA step, elsewhere prints a mmol-scale amount;
- protocol if the "1 mmol/L" recurs;
- OPD-limited if GO / MGO / DA dip when the osones peak.

**Test:**
1. Hofmann 2000 (same lab, GO / MGO SIDA with [13C4]butane-2,3-dione) prints the identical "(2 mL; 1 mmol/L)", and
   "1 mol/L" for the osone HPLC step, as Hofmann 1999 does. The 1 mmol/L is therefore the lab's written protocol,
   copied between papers. It may still be a copied misprint, but it is not a one-off.
2. Shape: 3-deoxypentosone peaks at 10 min (BPA 2304 -> 2438 -> 1192). MGO rises monotonically through it
   (0.345 -> 0.582 -> 0.763 umol/mmol at 5 / 10 / 20 min). DA is flat (BPA 35 / 34 / 25 / 16). Only GO falls (133 -> 18),
   and GO's fall is also chemically expected.
3. Budget: a quarter-pot holds about 13 umol of 3-DP at 10 min (BPA 2438 / CD 256 x 6 mmol / 4) plus 1-2.5 umol MGO,
   against 2 umol OPD. If the amount were real, OPD would be heavily limited. Because the IS ([13C4]-DA) competes for
   the same reagent, the bias would partly cancel for DA and not for MGO.

**Verdict:** OPEN. P(misprint for 1 mol/L) 0.5. The absence of a dip is weak evidence against limitation (P 0.6 that
it was not limiting).
- Decisive: the recovery of [13C4]-DA (absolute IS peak area vs the calibration mixtures), which only the authors
  have. Or the method reference [4] / [14] (Hofmann 1998, Hofmann et al. 1999), if either prints a recovery test with
  osones present.
- The question also covers hofmann2000 Fig. 3 GO / MGO (the same text).

```yaml
  - id: CF-067
    status: open
    explanations:
      - "Figure vs table: the figure was drawn from a different analysis or replicate set (or carries an unstated species-specific correction): P 0.6."
      - "(rejected 2026-11-01) A slip in Table 3's glycolaldehyde CD factor: the figure / table ratio is a constant 0.88 +- 0.02 at 10 / 60 / 120 min, not a power of 2; MGO on the same axis agrees; Rf 0.96 cannot produce it."
      - "OPD: '2 mL; 1 mmol/L' is the lab's written SIDA protocol (identical in hofmann2000), either a copied misprint for 1 mol/L (P 0.5) or real (then OPD-limited, and the bias does not cancel for MGO against a [13C4]-DA IS)."
    tests:
      - "2026-11-01: glycolaldehyde is quantified as its O-ethyloxime (ethoxamine, [13C4]-3-hydroxy-2-butanone IS), not with OPD: the OPD question does not touch the glycolaldehyde rows. MGO rises monotonically through the 10-min 3-deoxypentosone peak (no competition dip): weak evidence against OPD limitation."
    remedy: "Unchanged for glycolaldehyde (table value, 0.066 dex). MGO / GO / DA of Hofmann 1999 Fig. 10 and Hofmann 2000 Fig. 3: carry 'opd_amount_uncertain' (lower-bound alternative) until the IS recovery or the method reference decides."
    affects: [hofmann1999_fig10_xylose_ala_c2_c4_fragments, hofmann2000_fig3_dicarbonyls_glucose_phe]
```

## 8. CF-070: phenylalanine degraded in hofmann2000 Table 3 (MGO + Phe) vs hofmann2000b Table 3 (glucose + Phe)

**Rule:** the common rule, at hofmann2000's printed precision (integer %). For a common chemical step, equal Phe loss
must also be chemically possible for both donors.

**Test (cells from both PDFs):**

| t (min) | 2000: degraded % | 2000b: 1000 - unreacted (umol) | 2000b rounded | same digit? |
|---|---|---|---|---|
| 10 | 9 | 80 (920) | 8 | no (1 point) |
| 30 | 18 | 175 (825) | 17.5 -> 18 | yes |
| 60 | 36 | 360 (640) | 36 | yes |
| 120 | 49 | 490 (510) | 49 | yes |
| 240 | 58 | -- | -- | -- |
| 300 | 61 | 602 (398) | 60 | no (1 point) |
| 500 | 69 | -- | -- | -- |

- 3 of 5 cells match to the printed digit, so the rule outcome is DUPLICATE.
- Coincidence check: independent runs of one recipe scatter 7-16 % (CF-075). P(|diff| <= 0.5 point) is about 0.08 at
  36 % and about 0.06 at 49 %, so all five within 1 point gives P ~ 1e-3 for two runs of ONE recipe.
- Chemistry test: these are two different donors. Methylglyoxal (100 mM) forms imines with an alpha-amino acid orders
  of magnitude faster than glucose (open-chain fraction ~0.002 %), yet the two pots lose Phe identically from 10 min on.
- Each pot's own products close differently: glucose pot ARP-Phe 42 umol of 80 lost at 10 min; MGO pot PA + PAA
  16.3 umol of 90.

So a common step is implausible: P(coincidence + common step) ~ 0.1.

**Owner:**
- For the glucose pot: 2000b prints the assay (amino-acid analyser, ninhydrin), which reads Phe and ARP-Phe in the
  same chromatogram. Its 3-figure umol column feeds its yield column, which reproduces (52.5 % = 42 / 80). Hofmann 2000
  "never says how phenylalanine was measured".
- Against: hofmann2000 prints 7 times (adds 240 and 500 min). So it is not a subset of 2000b; both are subsets of one
  series.
- P(owner = glucose pot, 2000b) 0.7.

**Verdict:** EXPLAINED -> duplicate_confirmed. P 0.9 that it is one series; P 0.7 that it belongs to 2000b.

```yaml
  - id: CF-070
    status: duplicate_confirmed
    sources: [hofmann2000b, hofmann2000]      # 2000b first: the probable owner (assay printed; same chromatogram as ARP-Phe)
    tests:
      - "2026-11-01, both PDFs: 3 / 5 shared times identical to hofmann2000's printed digit (18 / 36 / 49 %), 2 within 1 point (9 vs 8.0, 61 vs 60.2); P ~1e-3 for two runs of one recipe; equal loss from MGO and glucose is chemically implausible (MGO far faster at forming the imine)."
    remedy: >-
      copy_of: hofmann2000_table3_mgo_phe_timecourse (PHE rows) -> hofmann2000b_table3_arp_phe_formation (PHE rows);
      the 240 / 500-min cells of hofmann2000 are extra points of the same series (glucose pot, P 0.7), never MGO-pot data.
      The MGO pot has no phenylalanine measurement. The PA / PAA rows of hofmann2000 Table 3 are unaffected.
```

## 9. CF-074: PA / PAA at 120 min, hofmann2000 Fig. 2 vs hofmann2000b Table 2 expt II (Glc / Phe, air)

**Rule:** the common rule (needs >= 3 cells). Only 2 exist (PA 1.399 vs 1.4; PAA 1.797 vs 1.8), so the rule outcome is
UNDECIDED.

**Likelihood check:**
- The Fig. 2 points sit on a 0.1 grid (the figure was "likely plotted from values rounded to 0.1"), so a match means
  landing in the same 0.1 bin.
- With a between-run SD of ~10-15 % (~0.15-0.25 umol/mmol for the difference), P(same bin) is about 0.2 per species,
  or about 0.05-0.08 for both.
- Same paper, other times: Fig. 2 at 30 min differs from Table 2 of hofmann2000 by 18 % (PA). 2000b Fig. 1 at 120 min
  (1.21 / 1.52) differs from Table 2 by 7-16 % (CF-075). Runs of this recipe in this lab do scatter.
- With a prior of 0.2-0.3 for reuse: **P(one run printed twice) ~ 0.8**.

**Control (pre-stated in the triage):**
- 2000b's ARP-Phe air / argon PA ratio is 5.5 / 0.6 = 9.2 (PAA 3.0 / 0.3 = 10), against Glc / Phe 1.4 / 0.4 = 3.5.
- log10(9.2 / 3.5) = 0.42 dex > 0.15, so by the criterion B29's 3.5 is flagged as possibly comparing a reflux pot
  (air) with a closed argon vial.
- Caveat: the control is not clean, since B29's premise is that ARP and sugar differ in O2 sensitivity.
- Same-lab, same-vessel reference: hofmann2000 Table 4 (closed vial, 60 min, dicarbonyl + Phe). Argon / air PA is
  11.1 / 10.2, 9.2 / 8.3, 8.4 / 7.8 and 8.1 / 6.9, so PA is insensitive to O2 with a dicarbonyl donor. The O2 effect on
  PA therefore runs through the Amadori / sugar step, which is what B29 models. That does not remove the vessel
  confound on the 3.5.

**Verdict:** EXPLAINED (one run, P 0.8). The B29 flag stays.

```yaml
  - id: CF-074
    status: explained
    tests:
      - "2026-11-01: only 2 cells (rule: undecided); Fig. 2 is plotted on a 0.1 grid, so P(both in the same bin by chance) ~0.05-0.08 given the lab's 7-18 % run scatter (Fig. 2 vs Table 2 of hofmann2000 at 30 min; 2000b Fig. 1 vs Table 2): P(one run) ~0.8."
      - "Control: ARP-Phe air / argon PA 9.2 (PAA 10) vs Glc / Phe 3.5: 0.42 dex apart (> 0.15), so the Glc / Phe ratio is kept only with the vessel flag. hofmann2000 Table 4 (closed vial): argon / air PA 1.09-1.17 with dicarbonyl donors."
    remedy: "Unchanged (count once; B29's 3.5 flagged 'possible_vessel_comparison'); B29's ARP-Phe 9.2 is the same-vessel pair and needs no flag."
```

---

## 10. Systematic scan: same-lab "value + SD" fingerprints

**Method.**
- `scratchpad/c1/fingerprint_scan.py` extracts every printed `value +- SD` pair (+-, ±, +/-) from all extraction
  dossiers and keyed fit tables (`data/lit/fit_rows/**`, including the `sd:` fields), normalised as printed (no unit
  conversion).
- It keeps pairs with >= 4 significant digits, drops fingerprints shared by > 6 documents (quoted constants), and lists
  document pairs sharing >= 3 (and, separately, exactly 2).
- Lab per dossier comes from the ledgers.
- 12,209 distinct fingerprints; 53 pairs with >= 3 shared.
- The measurement ledgers record WHERE the data are, never the values, so the values are scanned in the dossiers they
  point to.
- A second, value-only scan (`value_only_scan.py`, numbers with >= 4 significant digits shared by two dossiers of one
  lab) covers papers that print no SD (Leahy, Hofmann). Output: `value_only_hits.txt`.

**Hits that are reused data series (excluding the cluster 1 pairs above):**

| # | pair | evidence | already known? | P (reuse) | proposed action |
|---|---|---|---|---|---|
| H1 | bornhorst2017 (Part I) -> bornhorst2017b (Part II) 90 C rows | 17 value +- SD cells identical (68.2 +- 1.3, 39.5 +- 13.1, 28.3 +- 5.3, ...); one SE differs (L* k 6.0 +- 4.3 vs 6.0 +- 4.6) | ledger note ("90 C rows imported from Part I"); **not in the register** | 0.98 | Register as a duplicate (copy_of bornhorst2017b 90 C rows -> bornhorst2017_tab1); importance 0 (no wave reads either) |
| H2 | devleeschouwer2009 (Part II) Table 3 footnote-c cells -> DV 2008 / 2009b | 3.57, 0.10, 1.70, 1.23, 3.90, 26.43, 117.5, 105.7, 180.3, 105.4, 113.2, -6.7 and moisture 14.53 +- 0.08 | Part II dossier says "a second printing"; not in CF-040 | 0.99 | Fold into CF-040 (done above) |
| H3 | FUR / NF rows of wang2026b (4), zhou2022 (3-4) and hofmann1998 T5 (2) keyed in BOTH deoxyosone_submodel.yml and thiol_submodel.yml | identical value + sd, same ledger ids | handled: `fit_joint_deoxyosone_thiol.py::dedupe` counts them once (design review B7) | < 0.05 that some loader double-counts | none: the only other loader of both tables (`scripts/diagnostics/deoxypentosone_partition.py`) reads two named pots, not these rows |
| H4 | kang2026_SI <-> zhai2023foodchem (76 cells; cross-lab) | Kang's SI S4 reprints Zhai 2023 Table 1 | known (ledger original_of_kang2026_si_table_s4; CF-002 / 055 / 056) | 0.99 | none |
| H5 | zhang1988 / zhang1989 -> zhang1991b Table II totals (3782.6, 815.8, 2879.6, 2561.1) | value-only scan | known (registered, sulfur ledger duplicate_of) | 0.99 | none |
| H6 | thiol_hofmann_t9.yml <-> thiol_submodel.yml (MP rows) | identical rows | known: t9 status `merged` | -- | none |

**Inspected and rejected** (dossier prose quoting another paper's constant, physical constants, molecular weights, DOIs,
temperatures "120 +- 1 C", or chance matches across different compounds):
- knol2005 / knol2010, devleeschouwer / knol, claeys2005 / knol2005, martins2005 / piornos2025 / yu2018 / yu2020: the
  dossiers quote Martins' or De Vleeschouwer's constants as comparators;
- chan1994b / cremer2000, chan1994 / chan1994b (molar loadings);
- feng2022 / zhai2023foodchem and zhou2022 / zhou2023: comparison tables in the dossiers;
- bi2020 / bi2022: bi2022 cites Bi 2020's raw-flour levels, which no matrix entry of bi2022 uses as its own data (P
  reuse 0.1);
- kocadagli2016foodchem / jafc (a quoted ratio);
- deng2021 / deng2022, deng2022 / zhao2017, raza2020 / wang2026, sen2022 / xu2024, du2020 / zhao2017 / fan2018,
  Zhang2024 / zhang2024b, xu2007 / xu2010, lau2003 / pandit2006, hidalgo2016 / zamora2015 (not a duplicate, sec. 4).

**What the scan cannot see:**
- reprints with no printed SD and fewer than 4 significant digits (CF-074-type single cells);
- figure-only reprints (martins2003 / 2003b);
- re-plots in different units.

The value-only scan catches the Leahy-type case (15 shared numbers). The repo-side finding (one PDF dossiered under
several stems, CF-018) suggests a cheap gate: refuse two dossiers that cite the same `data/articles/*.pdf` or the same
DOI unless one declares `same_pdf_as`.

**Bottom line of the scan:** no unregistered same-lab reprint of a FITTED series was found. The only unregistered hit
is H1 (Bornhorst, unused). Single-lab mechanism evidence is not thinner than counted beyond what cluster 1 already
records (Leahy pH 9 arm; the Hofmann Phe series and Glc / Phe 120-min cell).
